import numpy as np
from astropy.io import fits
import os, sys
from os.path import dirname, abspath
from readsnapsgl import readsnap
import argparse
import yaml
import re
d = dirname(dirname(abspath(__file__)))
sys.path.append(d)
import pymgal



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



# The function which handles the data and creates the projections
def project_to_fits(object_dir, coords, config):
    handle_input(config)
   
    # Initialize some variables
    out_val = config["out_val"].lower()  
    sim_file_path = config["sim_file"]
    

    # Prepare SSP models, filters, and dust function
    sspmod = pymgal.SSP_models(config["SSP_model"], IMF=config["IMF"], has_masses=True)  # Prepare SSP models
    filters = pymgal.filters(f_name=config["filters"]) #eg: 'sloan_r', 'sloan_u', 'sloan_g', 'sloan_i', 'sloan_z', 'wfc3_f225w', 'wfc3_f606w', 'wfc3_f814w'])  # Load the filters
    dustf = None if config["dust_func"] is None else getattr(pymgal.dusts, config["dust_func"].lower())()   #Dust function
    print("Dust function being used: ", dustf) 
            
    snapname = os.path.basename(sim_file_path).replace('.hdf5', '') if ".hdf5" in os.path.basename(sim_file_path) else os.path.basename(sim_file_path)

    # Convert any custom projection vector to a numpy array and combine with the random arrays
    chosen_proj_vecs=[]
    if config["proj_vecs"] is not None:
        chosen_proj_vecs = [np.array(elem) if isinstance(elem, list) else elem for elem in config["proj_vecs"]]
    
    # Set up projection angles
    random_proj_vecs = random_projections(config["num_random_proj"])# random_proj_dict.get(snapname)
    projections = chosen_proj_vecs + random_proj_vecs + config["proj_angles"]
    cc, rr = coords[:3], coords[-1]
    
    head = readsnap(sim_file_path, 'HEAD')
    
    # Load the data, handle some user input, and calculate energy
    simd = pymgal.load_data(sim_file_path, snapshot=True, center=cc, radius=rr)
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


        # If the projection direction has type np.ndarray, change the format; otherwise DS9 can't open it
        if (isinstance(proj_direc, np.ndarray)):
            proj_direc = f"({';'.join(map(str, proj_direc))})"
        if (isinstance(proj_direc, list)):
            proj_direc = f"({'°,'.join(map(str, proj_direc))}°)"

        prefix = config["prefix"]
        prefix = prefix + "-" if prefix not in {None, ''}  else ""
        pj.write_fits_image(os.path.join(config["output_dir"], f"{prefix}{snapname}-{proj_direc}.fits"), overwrite=True)



if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="A python script to create projections of galaxies from hydrodynamical simulations.")
    
    parser.add_argument('sim_file', type=str, help="The simulation snapshot file you're using for your projection. Make sure to include the path as well as the filename.")
    parser.add_argument('config_file', type=str, help="The config.yaml file you're using to define your parameters. Make sure to include the path as well as the filename.")
    parser.add_argument('output_dir', type=str, help="Directory path where you want your files to be written.")
    parser.add_argument('--coords', type=float, nargs=4, help="The coordinates [x, y, z, r] of the projection region. (x, y, z) is the centre of the projection region and r is its radius")
    parser.add_argument('--prefix', type=str, help="If you want your output file to have a prefix (e.g. if you want to name this object belonging to these snaps. If null, no prefix will be included.")
    parser.add_argument('--filters', type=str, nargs='*', help="The filters you want for your mock observations")
    parser.add_argument('--num_random_proj', type=int, help="The number of random projections you want. Can be 0. NOTE: setting num_random_proj = 0, proj_vecs = null, and proj_angles = null will cause an error.")
    parser.add_argument('--proj_vecs', type=str, nargs='*', help="A list projection vectors. Can be 'x', 'y', 'z' (principal axes), and/or a custom [x, y, z] in cartesian. Can be null.")
    parser.add_argument('--proj_angles', type=float, nargs='*', help="A list of angles [alpha, beta, gamma] (degrees) used to rotate around the x, y, z axes respectively. Does the same task as proj_vecs, so you can use either/both. Can be null.")
    parser.add_argument('--out_val', type=str, help="Whether you want the output values to be 'flux', 'magnitude', or 'luminosity'")
    parser.add_argument('--mag_type', type=str, help="If out_val = 'magnitude', choose 'AB', 'vega', 'solar', or 'apparent' magnitude. If out_val != 'magnitude', this does nothing.")
    parser.add_argument('--rest_frame', type=lambda x: (x.lower() == 'true'), help="If True, results are in rest frame. Otherwise, they are in observer's frame.")
    parser.add_argument('--AR', type=float, help="Angular resolution (arcsec). If null, it is automatically calculated. If AR is null and npx is 'auto', npx defaults to 512")
    parser.add_argument('--SSP_model', type=str, help="The simple stellar population (SSP) model you want to use. Full list can be found at pymgal/models")
    parser.add_argument('--IMF', type=str, help="The initial mass function (IMF) assumed by the SSP model. Can be either 'chab', 'krou', or 'salp' (ie Chabrier, Kroupa, or Salpeter)")
    parser.add_argument('--dust_func', type=str, help="The assumed dust function. Can be either null, 'charlot_fall' or 'calzetti' (ie Charlot and Fall (2000) or Calzetti et al (2000))")
    parser.add_argument('--ksmooth', type=float, help="A smoothing parameter of at least 0 used in kNN Gaussian smoothing. The greater ksmooth, the smoother.")
    parser.add_argument('--z_obs', type=float, help="The distance of the observation from the observer's POV. This DOES NOT affect age/evolution, only the apparent distance. If null, it defaults to max(0.05, sim z)")
    parser.add_argument('--npx', type=int, help="The number of pixels in the output image. Can also be set to 'auto' which will automatically decide the pixel number to include all particles")
    parser.add_argument('--zthick', type=float, help="The thickness cut ([center-zthick, center+zthick]) in the projection direction in kpc/h. If null, all the data is used (ie no cut is applied)")
    parser.add_argument('--outmas', type=lambda x: (x.lower() == 'true'), help="Output the mass map corresponding to your data")
    parser.add_argument('--outage', type=lambda x: (x.lower() == 'true'), help="Output the age map corresponding to your data")
    parser.add_argument('--outmet', type=lambda x: (x.lower() == 'true'), help="Output the metallicity map corresponding to your data")
    args = parser.parse_args()

    cwd = os.getcwd()
    config_file = os.path.join(cwd, "config.yaml")
    config=None  # Initialize
    # Exit if the file doesn't exist or is unreadable. Otherwise load it.
    if not os.path.isfile(config_file) or not os.access(config_file, os.R_OK):
        sys.exit(f"Error: Cannot read file {config_file}" if os.path.isfile(config_file) else f"Error: File {config_file} does not exist.")
    else:
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
    config = merge_settings(config, args)
    project_to_fits(config["sim_file"], config["coords"], config)

