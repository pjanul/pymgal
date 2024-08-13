import yaml
import subprocess
import numpy as np
from astropy.io import fits
import os, sys
import argparse
import yaml
import re
global start_index, end_index, snaps, code, cat_dir, sim_dir, output_dir, config_file_path

##################################################################################################################################################################################################################

################################################## This script helps you generate projections from the 300 project simulation of galaxy clusters. ###########################################################


# If you're not a member of the 300 project but have AHF catalogue data, you might still be able to modify this script to get the positions of halos and their progenitors
# Below you'll find some parameters including catalogue directory, simulation directory, etc.
# Modify these to suit your needs and then run the script at the command line.
# We assume the following directory structure.

#       ├── sim_dir
#            ├── NewMDCLUSTER_0001
#                ├── snap_128.hdf5
#                ├── snap_127.hdf5
#                ├── more snapshot files
#            ├── more NewMDCLUSTER_XXXX folders containing snapshot files
#
#        ├── cat_dir
#            ├── GIZMO_R200c_snaps_128-center-cluster.txt
#            ├── NewMDCLUSTER_0001
#                ├── GIZMO-NewMDCLUSTER_0001.snap_128.z0.000.AHF_halos
#                ├── GIZMO-NewMDCLUSTER_0001.snap_128.z0.000.AHF_mtree_idx
#                ├── AHF_halos, and AHF_mtree_idx files for other snap numbers...
#            ├── more NewMDCLUSTER_XXXX folders containing .AHF_halos and .AHF_mtree_idx files


# IMPORTANT NOTE: You should move this script inside the inner PyMGal directory and run it from there. It calls "python ./generate_mocks.py", so it needs to be in the same directory as generate_mocks.py.
# After moving the script to the proper location, enter  "python the300_helper.py <config_filepath> <output_dir>" into the command line.

#################################################################################################### Here are the parameters you should modify #####################################################################################################

# Example input values
start_index = 300                                                                           # Index of the first cluster you want to project
end_index = None                                                                            # Index of the last cluster you want to project. If None, only start_index will be projected
snaps = ["snap_128", "snap_127"]
code = "GIZMO"                                                                              # What is your simulation code? Can be GadgetX, GIZMO, GIZMO_7k, etc
cat_dir = "/niftydata/TheThreeHundred/data/catalogues/AHF/GIZMO"                            # The directory to your AHF catalogue data. You'll need a one single center-clusters.txt file and .AHF_halos, .AHF_mtree_idx for each cluster/redshift pair
sim_dir = "/GIZMO-SIMBA"                                                                    # The directory to your simulation snapshot data

#################################################################################################  Below is code you probably don't need to touch ########################################################################################################



def generate_object_names(start_index, end_index):
    # Generate cluster names given the first and last index
    object_names = []
    if end_index is None:
        return [f"NewMDCLUSTER_{start_index:04d}"]
    for i in range(start_index, end_index + 1):
        object_name = f"NewMDCLUSTER_{i:04d}"
        object_names.append(object_name)
    
    return object_names


def get_progenitors(snaps, cluster_name, sim_dir, cat_dir, code, config, output_dir, sim_file_ext):

    snapnum_list = [int(snap.replace("snap_", "")) for snap in snaps]

    # Load group information
    groupinfo = np.loadtxt(os.path.join(cat_dir, f"{code}_R200c_snaps_128-center-cluster.txt"))

    pregen = np.zeros((groupinfo.shape[0], 129), dtype=np.int64) - 1
    pregen[:, 128] = np.int64(groupinfo[:, 1])

    try:
        cluster_num = cluster_name.split('_')[1]
    except:
        raise ValueError("If you're using AHF halos and merger trees to get coordinates, they must end in _X, where X is some sequence of digits (e.g. NewMDCLUSTER_0324)")

    # See how far back your progenitors need to go
    earliest_snap_num = min([int(snapname.split('_')[-1]) for snapname in snaps])

    lp = [int(cluster_num) - 1]

    # Get progenitors.
    for i in np.arange(128, earliest_snap_num, -1):
        exts = '000' + str(i)
        exts = exts[-3:]  # Ensures it's 3 digits

        # Construct the path to the merger tree index file directory
        mtree_dir = os.path.join(cat_dir, cluster_name)

        # List all files in the directory and filter for the correct pattern
        files = os.listdir(mtree_dir)
        mtree_files = [f for f in files if f.startswith(f'{code}-{cluster_name}.snap_{exts}.z') and f.endswith('.AHF_mtree_idx')]

        if mtree_files:
            mtree_index_file = os.path.join(mtree_dir, mtree_files[0])  # Choose the first match if multiple files exist
            mtid = np.loadtxt(mtree_index_file)
            print("AHF merger tree index file successfully found: ", mtree_index_file)
        else:
            raise ValueError(f'Cannot find merger tree index file for snap_{exts} in {mtree_dir}')

        pregen[lp, i - 1] = mtid[mtid[:, 0] == pregen[lp, i], 1]

    return pregen

def get_coords(snapname, cluster_name, sim_path, cat_path, Code, pregen, sim_file_ext):

    cluster_num = cluster_name[-3:]

    lp = [int(cluster_num) -1]
    sn = np.int64(snapname[-3:])
    hid= pregen[lp,sn]

    halo_dir = os.path.join(cat_path, cluster_name)
    files = os.listdir(halo_dir)
    halo_files = [f for f in files if f.startswith(f'{Code}-{cluster_name}.{snapname}.z') and f.endswith('.AHF_halos')]



    if halo_files:
        halo_file = os.path.join(halo_dir, halo_files[0])  # Choose the first match if multiple files exist
        halo = np.loadtxt(halo_file)
        print("AHF halo file successfully found: ", halo_file)
    else:
        raise ValueError('Cannot find halo file  %s' % halo_file)

    shid = np.where(np.int64(halo[:,0]) == hid)[0]
    if len(shid) == 0:
        raise ValueError('Cannot find halos!! %s' % snapname)

    for j in shid:
        cc=halo[j,5:8]; rr = halo[j,11]
        coords = list(np.append(cc, rr))

    return coords



def process_snaps(snaps, cluster_name, sim_dir, cat_dir, code, config, output_dir, config_file_path):

    sim_file_ext = '' if "gadget" in code.lower() else '.hdf5'

    output_dir = os.path.join(output_dir, cluster_name)
    pregen = get_progenitors(snaps, cluster_name, sim_dir, cat_dir, code, config, output_dir, sim_file_ext)
    for snap in snaps:
        sim_file = os.path.join(sim_dir, cluster_name, snap + sim_file_ext)
        coords = get_coords(snap, cluster_name, sim_dir, cat_dir, code, pregen, sim_file_ext)
        coords_str = coords = [str(coord) for coord in coords]   
        
        # Run generate_mocks.py at the command line
        command = ["python", "generate_mocks.py", sim_file, config_file_path, output_dir]
        command.extend(["--coords"] + coords_str)
        command.extend(["--prefix", cluster_name])
        subprocess.run(command)


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Generate projections from the 300 project simulation of galaxy clusters.")
    parser.add_argument('config_dir', type=str, help="Path to the config file")
    parser.add_argument('output_dir', type=str, help="Path to the output directory")

    args = parser.parse_args()

    config_file_path = args.config_dir
    output_dir = args.output_dir

    # Ensure the output directory exists
    if not os.path.isdir(output_dir):
        raise ValueError("The output directory does not exist")

    # Read the YAML config file
    with open(config_file_path, 'r') as file:
        config = yaml.safe_load(file)

    # Generate object names
    object_names = generate_object_names(start_index, end_index)

    # Iterate through each object name and process the snaps
    for cluster_name in object_names:
        process_snaps(snaps, cluster_name, sim_dir, cat_dir, code, config, output_dir, config_file_path)

if __name__ == "__main__":
    main()
