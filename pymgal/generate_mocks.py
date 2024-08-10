import numpy as np
from astropy.io import fits
import os, sys
from os.path import dirname, abspath
from readsnapsgl import readsnapgd, readhdf5head, readsnap
import argparse
import yaml
import re
import json
import multiprocessing
d = dirname(dirname(abspath(__file__)))
sys.path.append(d)
import pymgal


# Worker functions for parallel processing
def process_object(object_dir, config):
    print(f"-------------------- {object_dir} ---------------------")
    project_to_fits(object_dir, config)

def process_snapshot(params):
    object_dir, snapname, config = params
    config["snaps"] = [snapname]  # Process one snapshot at a time
    print(f"-------------------- {object_dir} ---------------------")
    project_to_fits(object_dir, config)



# Generate a list of cluster numbers such as "0001", "0002", etc.
def create_number_strings(start, end):
    number_strings = []
    if end is None:
        number_strings = [str(start).zfill(4)]
    elif end < start:
        raise ValueError("The starting cluster cannot have a greater index than the end cluster.")
    else:
        for num in range(start, end + 1):
            num_string = str(num).zfill(4)
            number_strings.append(num_string)

    return number_strings


# Generate n random projection vectors
def random_projections(n_samples):
    points = []
    for _ in range(n_samples):
        # Generate random point on unit sphere using spherical coordinates
        theta = np.random.uniform(0, 2*np.pi)
        phi = np.random.uniform(0, np.pi)

        x = np.sin(phi) * np.cos(theta)
        y = np.sin(phi) * np.sin(theta)
        z = np.cos(phi)

        points.append([round(x, 2), round(y, 2), round(z, 2)])
    points = [np.array(element) for element in points]
    return points

# Overrule the config.yaml file if the user inputs something at the command line
def merge_settings(config, args):
    # Update settings from config with values from args
    for key, value in vars(args).items():
        if value is not None:
            config[key] = value
    return config


def handle_input(config):
    # Handle some invalid inputs that aren't handled elsewhere
    if config["out_val"].lower() not in {"flux", "magnitude", "luminosity"}:
        raise ValueError(f"Invalid out_val: {config['out_val']} Expected 'flux', 'magnitude', or 'luminosity'.") 

    if (config["mag_type"].lower() not in {"vega", "solar", "apparent", "ab"}) and (config["out_val"]  == "magnitude"):
        raise ValueError(f"Invalid mag_type: {config['mag_type']}. Expected 'vega', 'solar', 'apparent', or 'AB'.")
    


def get_progenitors(snaps, object_dir, sim_path, cat_path, Code, sim_file_ext=''):
    
    snapnum_list = [int(snap.replace("snap_", "")) for snap in snaps]
    
    # Load group information
    groupinfo = np.loadtxt(os.path.join(cat_path, f"{Code}_R200c_snaps_128-center-cluster.txt"))
    pregen = np.zeros((groupinfo.shape[0], 129), dtype=np.int64) - 1
    pregen[:, 128] = np.int64(groupinfo[:, 1])

    try:
        cluster_num = object_dir.split('_')[1]
    except:
        raise ValueError("If you're using AHF halos and merger trees to get coordinates, they must end in _X, where X is some sequence of digits (e.g. NewMDCLUSTER_0324)")

    # See how far back your progrenitors need to go
    earliest_snap_num = min([int(snapname.split('_')[-1]) for snapname in snaps])

    lp = [int(cluster_num) -1]

    # Get progenitors.
    for i in np.arange(128, earliest_snap_num, -1):
        exts='000'+str(i)

        head_path = os.path.join(sim_path, object_dir, 'snap_' + exts[-3:] + sim_file_ext)
        head = readsnap(head_path, 'HEAD') 

        if head.Redshift<0:
            head.Redshift = 0.0000
        mtree_index_file = os.path.join(cat_path, object_dir, Code + '-' + object_dir + '.snap_' + exts[-3:] + '.z' + ("{:.3f}".format(head.Redshift, 9))[:5] + '.AHF_mtree_idx')

        if os.path.isfile(mtree_index_file):
            mtid=np.loadtxt(mtree_index_file)
            print("AHF merger tree index file successfully found: ", mtree_index_file)
        else:
            raise ValueError('Cannot find merger tree index file %s' % mtree_index_file)


        pregen[lp,i-1]=mtid[mtid[:,0]==pregen[lp,i],1]
    return pregen


def get_coords(snapname, object_dir, sim_path, cat_path, Code, pregen, custom=False, sim_file_ext=''):

    if custom:

        print("Using user-defined coordinates")
        coords_path = os.path.join(cat_path, object_dir, object_dir + "_" + snapname + "_coords.npy")
        arr = np.load(coords_path, allow_pickle=True)
        cc, rr = arr[:3], arr[-1]
        return cc, rr

    cluster_num = object_dir[-3:]

    lp = [int(cluster_num) -1]
    sn = np.int64(snapname[-3:])
    hid= pregen[lp,sn]

    head_path = os.path.join(sim_path, object_dir, snapname + sim_file_ext)
    head = readsnap(head_path, 'HEAD')
    if head.Redshift<0:
        head.Redshift = 0.0000
    halo_file=os.path.join(cat_path, object_dir, Code + '-'+ object_dir +'.'+snapname+'.z'+("{:.3f}".format(head.Redshift,9))[:5]+'.AHF_halos')

    if os.path.isfile(halo_file):
        halo=np.loadtxt(halo_file)
        print("AHF halo file successfully found: ", halo_file)
    else:
        raise ValueError('Cannot find halo file  %s' % halo_file)

    shid = np.where(np.int64(halo[:,0]) == hid)[0]
    if len(shid) == 0:
        raise ValueError('Cannot find halos!! %s' % snapname)

    for j in shid:
        cc=halo[j,5:8]; rr = halo[j,11]
        print("Centre, radius", cc, rr)

    return cc, rr



# The function which handles the data and creates the projections
def project_to_fits(object_dir, config):

    handle_input(config)

    # Initialize some variables
    snaps = config["snaps"]
    Code = config["code"]
    out_val = config["out_val"].lower()  
    cat_path = config["cat_dir"] 
    sim_path = config["sim_dir"] 
    
    # Prepare SSP models, filters, and dust function
    sspmod = pymgal.SSP_models(config["SSP_model"], IMF=config["IMF"], has_masses=True)  # Prepare SSP models
    filters = pymgal.filters(f_name=config["filters"]) #eg: 'sloan_r', 'sloan_u', 'sloan_g', 'sloan_i', 'sloan_z', 'wfc3_f225w', 'wfc3_f606w', 'wfc3_f814w'])  # Load the filters
    dustf = None if config["dust_func"] is None else getattr(pymgal.dusts, config["dust_func"].lower())()   #Dust function
    print("Dust function being used: ", dustf) 
            
    sim_file_ext = ''
    if not os.path.isfile(os.path.join(sim_path, object_dir, 'snap_128')) and os.path.isfile(os.path.join(sim_path, object_dir, 'snap_128' + ".hdf5")):
        sim_file_ext = '.hdf5'
    elif os.path.isfile(os.path.join(sim_path, object_dir, 'snap_128')) and os.path.isfile(os.path.join(sim_path, object_dir, 'snap_128' + ".hdf5")):
        raise ValueError("Do not include multiple simulation snapshot files in the same directory.")


    pregen = None
    if not config["custom_coords"]:
        pregen = get_progenitors(snaps, object_dir, sim_path, cat_path, Code, sim_file_ext)
    
    for snapname in snaps:
        # Convert any custom projection vector to a numpy array and combine with the random arrays
        chosen_proj_vecs=[]
        if config["proj_vecs"] is not None:
            chosen_proj_vecs = [np.array(elem) if isinstance(elem, list) else elem for elem in config["proj_vecs"]]

        random_proj_vecs = random_projections(config["num_random_proj"])# random_proj_dict.get(snapname)
        projections = chosen_proj_vecs + random_proj_vecs + config["proj_angles"]
            
        cc, rr = get_coords(snapname, object_dir, sim_path, cat_path, Code, pregen, custom=config["custom_coords"], sim_file_ext=sim_file_ext) 

        head_path = os.path.join(sim_path, object_dir, snapname + sim_file_ext)
        head = readsnap(head_path, 'HEAD')
        
        # Load the data, handle some user input, and calculate energy
        simd = pymgal.load_data(head_path, snapshot=True, center=cc, radius=rr*1.4)
        mag_type = config["mag_type"].lower()
        is_vega, is_solar, is_apparent = (mag_type == "vega", mag_type == "solar", mag_type == "apparent")
        z_obs = config["z_obs"] if config["z_obs"] is not None else max(0.05, simd.redshift)
        mag = filters.calc_energy(sspmod, simd, dust_func=dustf, vega=is_vega, apparent=is_apparent, solar=is_solar, unit=out_val, rest_frame=config["rest_frame"])
                
        # Rotate and project
        for i, proj_direc in enumerate(projections):
            pj = None #Initialize
            if i==0:
                print("Projecting photons to %s" % proj_direc)
                pj = pymgal.projection(mag, simd, npx=config["npx"], unit=out_val, AR=config["AR"], redshift=z_obs, zthick=config["zthick"],
                                               axis=proj_direc, ksmooth=config["ksmooth"], outmas=config["outmas"], outage=config["outage"], outmet=config["outmet"])
                lsmooth = pj.lsmooth
                     
            # Avoid redundantly recomputing kNN distances by passing the pre-computed array 
            else:
                print("Projecting photons to %s" % proj_direc)
                pj = pymgal.projection(mag, simd, npx=config["npx"], unit=out_val, AR=config["AR"], redshift=config["z_obs"], zthick=config["zthick"],
                                           axis=proj_direc, ksmooth=config["ksmooth"], lsmooth=lsmooth, outmas=config["outmas"], outage=config["outage"], outmet=config["outmet"])


            # Create directory if it doesn't already exist
            output_dir = os.path.join(config["output_dir"], "maps", "CCD", Code, object_dir)
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)

            # If the projection direction has type np.ndarray, change the format; otherwise DS9 can't open it
            if (isinstance(proj_direc, np.ndarray)):
                proj_direc = f"({';'.join(map(str, proj_direc))})"
            if (isinstance(proj_direc, list)):
                proj_direc = f"({'°,'.join(map(str, proj_direc))}°)"
            pj.write_fits_image(os.path.join(output_dir, f"{object_dir}-{snapname}-{proj_direc}.fits"), overwrite=True)



if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="A python script to create projections of galaxy clusters using The 300 Project simulations.")
    parser.add_argument("output_dir", type=str, help="Directory path where you want your files to be written.")
    parser.add_argument("--config_file", type=str, default="./config.yaml", help="Path to your config file. If you have your own custom config file, enter it here. \
                        Note that some of the parameters in the config file can be overruled by passing arguments to this script.\n \
                        Default: ./config.yaml")
    parser.add_argument("--filters", nargs='+', help="A list of strings corresponding to the filters you'd like for your projections.")
    parser.add_argument("--snaps", nargs='+', help="A list of strings corresponding to the snap numbers you'd like for your projections. \n \
                        Note that snap_128 corresponds to redshift 0 and decreasing the snap number increases the redshift.")
    parser.add_argument("--proj_vecs", type=int, help="A list projection vectors. Can be 'x', 'y', 'z' (principal axes), or a custom [x, y, z] in cartesian. \n \
                        Can be null only if num_random_proj > 0.")
    parser.add_argument("--num_random_proj", type=int, help="Number of random projections that will be generated for each snap_num. \n \
                        Can be set to zero only proj_vecs is not null. If proj_vecs is not null and num_random_proj>0, all the vectors will be used" )
    args = parser.parse_args()


    config=None  # Initialize
    # Exit if the file doesn't exist or is unreadable. Otherwise load it.
    if not os.path.isfile(args.config_file) or not os.access(args.config_file, os.R_OK):
        sys.exit(f"Error: Cannot read file '{args.config_file}'." if os.path.isfile(args.config_file) else f"Error: File '{args.config_file}' does not exist.")
    else:
        with open(args.config_file, 'r') as f:
            config = yaml.safe_load(f)
    config = merge_settings(config, args)



    #cluster_nums = create_number_strings(config["start_index"], config["end_index"])
    # The list of objects you want to project. Directories must be named after the objects
    #for object_dir in config["obj_list"]:
    #    print(f"------------------- {object_dir} -------------------")
    #    project_to_fits(object_dir, config)

    # Determine parallelization strategy
    num_objects = len(config["obj_list"])
    num_snaps = len(config["snaps"])

    if num_objects > num_snaps:
        print("Pooling over objects. Using {config['ncpu']} CPUs")
        # Parallelize over objects
        with multiprocessing.Pool(config["ncpu"]) as pool:
            pool.map(lambda obj: process_object(obj, config), config["obj_list"])
    else:
        print(f"Pooling over snaps. Using {config['ncpu']} CPUs")
        # Parallelize over snapshots for each object
        params_list = [(obj, snap, config) for obj in config["obj_list"] for snap in config["snaps"]]
        with multiprocessing.Pool(config["ncpu"]) as pool:
            pool.map(process_snapshot, params_list)
