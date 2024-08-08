# README


## What is this repository for?

-   This package uses simple stellar synthesis models to generate observed galaxies from hydro-dynamical simulations.


## How do I get set up?

- Begin by downloading/cloning pymgal onto your own device and placing it in your directory of choice. After downloading, you will need to install the necessary dependencies and configure your database before you can run the code.

#### Dependencies 

To install the necessary dependencies, simply enter the outer pymgal directory and run the following at the command line.

- `pip install -r requirements.txt`


#### Database configuration


For a given cluster at a given redshift, you'll need the corresponding (1) snap_XYZ snapshot file (2) .AHF_halos file, (3) .AHF_mtree file, and (4) .AHF_mtree_idx files. You'll also need (5) a single center-cluster.txt file for your chosen simulation code. The snapshot files must be stored in the simulation directory (sim_dir), while the other files must be stored in the catalogue directory (cat_dir). If you're only interested in redshift z=0 (i.e. snap_128), you can omit the .AHF_mtree and .AHF_mtree_idx files.The contents of sim_dir and cat_dir must follow the format below. This example uses GIZMO data, but you can also use GadgetX if you'd like. 

<br />
```bash
├── sim_dir
    ├── NewMDCLUSTER_0001
        ├── snap_128.hdf5
        ├── snap_127.hdf5
        ├── more snap_ files 
    ├── More NewMDCLUSTER_ directories...

├── cat_dir 
    ├── GIZMO_R200c_snaps_128-center-cluster.txt
    ├── NewMDCLUSTER_0001
         ├── GIZMO-NewMDCLUSTER_0001.snap_128.z0.000.AHF_halos
         ├── GIZMO-NewMDCLUSTER_0001.snap_128.z0.000.AHF_mtree 
         ├── GIZMO-NewMDCLUSTER_0001.snap_128.z0.000.AHF_mtree_idx
         ├── More AHF_halos, AHF_mtree, and AHF_mtree_idx files...
    ├── More NewMDCLUSTER_ directories...
```  
<br />

#### The config.yaml file

- The config.yaml file contains modifiable parameters including filters, clusters, snapshots, and many more. To begin, make sure that sim_path and cat_path are correctly identified. You can play around with other parameters, but the default values should be enough to get you started. The file can be found at your/path/to/pymgal/pymgal (i.e. the inner pymgal directory). 


## How do I run the code?

Once everything is ready, you can begin generating mock observations. To do so, run the following at command line.

- **`python /your/path/to/pymgal/pymgal generate_mocks.py <output directory>`**


You can pass optional command line arguments to temporarily overrule the config file. You can even pass your own config file path as an argument and use that one instead.  To get more information about command line arguments, you can run the following. <br>

- **`python generate_mocks.py --help`**
<br />

## Who do I talk to?

-   Please report any issue to Weiguang Cui cuiweiguang@gmail.com.
-   Or report a bug through issues.

## Acknowledgement

-   This package borrowed a lot things from ezgal (<http://www.baryons.org/ezgal/>). Please make your acknowledgement to their work when you use this package.
