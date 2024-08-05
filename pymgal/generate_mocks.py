import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from astropy.io import fits
import os, sys
from os.path import dirname, abspath
from readsnapsgl import readsnapgd, readhdf5head
from PIL import Image
import argparse
import yaml
d = dirname(dirname(abspath(__file__)))
sys.path.append(d)
import pymgal


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

# The function which handles the data and creates the projections
def project_to_fits(cluster_num, config):

    # Initialize some variables
    snaps = config["snaps"]
    Code = config["code"]
    outiv = config["outiv"]  
    cat_path = config["cat_dir"] 
    sim_path = config["sim_dir"] 
    chosen_dustf = config["dust_func"]


    # Load group information
    groupinfo = np.loadtxt(os.path.join(cat_path, f"{Code}_R200c_snaps_128-center-cluster.txt"))
    pregen = np.zeros((groupinfo.shape[0], 129), dtype=np.int64) - 1
    pregen[:, 128] = np.int64(groupinfo[:, 1])

    halo = np.loadtxt(os.path.join(cat_path, f"NewMDCLUSTER_{cluster_num}", f"{Code}-NewMDCLUSTER_{cluster_num}.snap_128.z0.000.AHF_halos"))

    sspmod = pymgal.SSP_models(config["SSP_model"], IMF=config["IMF"], has_masses=True)  # Prepare SSP models
    filters = pymgal.filters(f_name=config["filters"]) #eg: 'sloan_r', 'sloan_u', 'sloan_g', 'sloan_i', 'sloan_z', 'wfc3_f225w', 'wfc3_f606w', 'wfc3_f814w'])  # Load the filters
    dustf = None if chosen_dustf is None else getattr(pymgal.dusts, chosen_dustf)()   #Dust function

    
    sim_file_ext = '.hdf5' if Code == 'GIZMO' else ''   # GadgetX simulation snapshot files have no file extension, while GIZMO has .hdf5

    for lp in [int(cluster_num)-1]:
        clnum='0000'+str(lp+1)
        clnum=clnum[-4:]
        cname = "NewMDCLUSTER_"+clnum+"/"

        # See how far back your progrenitors need to go
        earliest_snap_num = min([int(snapname.split('_')[-1]) for snapname in snaps])

        # Get progenitors.
        for i in np.arange(128, earliest_snap_num, -1):
            exts='000'+str(i)
            head_path = os.path.join(sim_path, f'NewMDCLUSTER_{cluster_num}', 'snap_' + exts[-3:] + sim_file_ext)
            head = readsnapgd(head_path, 'HEAD') if Code == 'GadgetX' else readhdf5head(head_path, 'HEAD') if Code == 'GIZMO' else None



            if head.Redshift<0:
                head.Redshift = 0.0000
            mtree_index_file = os.path.join(cat_path, cname + Code + '-' + cname[:-1] + '.snap_' + exts[-3:] + '.z' + ("{:.3f}".format(head.Redshift, 9))[:5] + '.AHF_mtree_idx')

            if os.path.isfile(mtree_index_file):
                mtid=np.loadtxt(mtree_index_file)
                print("AHF merger tree index file successfully found: ", mtree_index_file)
            else:
                raise ValueError('Cannot find merger tree index file %s' % mtree_index_file)


            pregen[lp,i-1]=mtid[mtid[:,0]==pregen[lp,i],1]



        for snapname in snaps:
            # Convert any custom projection vector to a numpy array and combine with the random arrays
            chosen_proj_vecs=[]
            if config["proj_vecs"] is not None:
                chosen_proj_vecs = [np.array(elem) if isinstance(elem, list) else elem for elem in config["proj_vecs"]]
            random_proj_vecs =  random_proj_dict.get(snapname)
            projections = chosen_proj_vecs + random_proj_vecs

            # Load halos and merger trees

            head_path = os.path.join(sim_path, f'NewMDCLUSTER_{cluster_num}', snapname + sim_file_ext)
            head = readsnapgd(head_path, 'HEAD') if Code == 'GadgetX' else readhdf5head(head_path, 'HEAD') if Code == 'GIZMO' else None


            if head.Redshift<0:
                head.Redshift = 0.0000
            sn = np.int64(snapname[-3:])
            hid= pregen[lp,sn]
            halo_file=os.path.join(cat_path, cname, Code + '-'+cname[:-1]+'.'+snapname+'.z'+("{:.3f}".format(head.Redshift,9))[:5]+'.AHF_halos')
            
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
                
                # Load simulation data
                simd = pymgal.load_data(head_path, snapshot=True, center=cc, radius=rr*1.4)

                # Calculate luminosity
                mag = filters.calc_energy(sspmod, simd, Ncpu=16, unit=outiv, rest_frame=True)

                # Rotate and project
                z_obs = config["z_obs"] if config["z_obs"] is not None else max(0.10, simd.redshift)

                # Rotate and project
                for i, proj_direc in enumerate(projections):

                    pj = None #Initialize
                    if i==0:
                        print("Projecting photons to %s" % proj_direc)
                        pj = pymgal.projection(mag, simd, npx=config["npx"], unit=outiv, AR=config["AR"], redshift=config["z_obs"], zthick=config["zthick"],
                                               axis=proj_direc, ksmooth=config["ksmooth"], outmas=config["outmas"], outage=config["outage"], outmet=config["outmet"])
                        lsmooth = pj.lsmooth
                     
                    # Avoid redundantly recomputing kNN distances by passing the pre-computed array 
                    else:
                        print("Projecting photons to %s" % proj_direc)
                        pj = pymgal.projection(mag, simd, npx=config["npx"], unit=outiv, AR=config["AR"], redshift=config["z_obs"], zthick=config["zthick"],
                                           axis=proj_direc, ksmooth=config["ksmooth"], lsmooth=lsmooth, outmas=config["outmas"], outage=config["outage"], outmet=config["outmet"])


                    # Create directory if it doesn't already exist
                    output_dir = config["output_dir"] + f"/maps/CCD/{Code}/NewMDCLUSTER_{cluster_num}"
                    if not os.path.exists(output_dir):
                        os.makedirs(output_dir)

                    # If the projection direction has type np.ndarray, change the format; otherwise DS9 can't open it
                    if (isinstance(proj_direc, np.ndarray)):
                        proj_direc = f"({';'.join(map(str, proj_direc))})"

                    pj.write_fits_image(output_dir + f"/NewMDCLUSTER_{cluster_num}-{snapname}-{proj_direc}.fits", overwrite=True)



if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="A python script to create projections of galaxy clusters using The 300 Project simulations.")
    parser.add_argument("output_dir", type=str, help="Directory path where you want your files to be written.")
    parser.add_argument("--config_file", type=str, default="./config.yaml", help="Path to your config file. If you have your own custom config file, enter it here. \
                        Note that some of the parameters in the config file can be overruled by passing arguments to this script.\n \
                        Default: ./config.yaml")
    parser.add_argument("--start_cluster", type=int, help="An integer between 1 and 324 which represents the index of the first cluster you'd like to project. \n  \
                        If you want to project a range of clusters in a single run, you must also specify the final cluster index via the --end_cluster argument.\n \
                        If you only want to project a single cluster, pass its index here and pass nothing to the --end_cluster argument.")
    parser.add_argument("--end_cluster", type=int, help="An integer between 1 and 324 which represents the index of the final cluster you'd like to project. \n \
                        If nothing is passed here, a single cluster specified by the --start_cluster argument will be projected. \n \
                        Default value: None.")
    parser.add_argument("--filters", nargs='+', help="A list of strings corresponding to the filters you'd like for your projections.")
    parser.add_argument("--snaps", nargs='+', help="A list of strings corresponding to the snap numbers you'd like for your projections. \n \
                        Note that snap_128 corresponds to redshift 0 and decreasing the snap number increases the redshift.")
    parser.add_argument("--proj_vecs", type=int, help="A list projection vectors. Can be 'x', 'y', 'z' (principal axes), or a custom [x, y, z] in cartesian. \n \
                        Can be null only if num_random_proj > 0.")
    parser.add_argument("--num_random_proj", type=int, help="Number of random projections that will be generated for each snap_num. \n \
                        Can be set to zero only proj_vecs is not null. If proj_vecs is not null and num_random_proj>0, all the vectors will be used" )
    args = parser.parse_args()

    # Overrule the config.yaml end_cluster if the user inputs a start cluster with no end
    if args.start_cluster and not args.end_cluster:
        args.end_cluster = None

    config=None  # Initialize
    # Exit if the file doesn't exist or is unreadable. Otherwise load it.
    if not os.path.isfile(args.config_file) or not os.access(args.config_file, os.R_OK):
        sys.exit(f"Error: Cannot read file '{args.config_file}'." if os.path.isfile(args.config_file) else f"Error: File '{args.config_file}' does not exist.")
    else:
        with open(args.config_file, 'r') as f:
            config = yaml.safe_load(f)
    config = merge_settings(config, args)


    # Create a random set of projection angles for each snap
    random_proj_dict = {}
    for snap_name in config["snaps"]:
        random_proj_dict[snap_name] = random_projections(config["num_random_proj"])

    # Define clusters and create the mock observations
    cluster_nums = create_number_strings(config["start_cluster"], config["end_cluster"])

    for cluster_num in cluster_nums:
        print(f"------------------- Cluster {cluster_num} -------------------")
        project_to_fits(cluster_num, config)
