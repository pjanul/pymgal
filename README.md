# README


## What is this repository for?

* PyMGal (pronounced py-em-gal) is a package that uses simple stellar synthesis models to generate observed galaxies from hydro-dynamical simulations.

## How do I get set up?

* Begin by downloading/cloning PyMGal onto your own device and placing it in your directory of choice. After downloading, you will need to install the necessary dependencies and configure your database before you can run the code.

#### Dependencies 

* To install the necessary dependencies, simply enter the outer pymgal directory and run the following at the command line.

    * `pip install -r requirements.txt`


#### The config.yaml file

* The config.yaml file contains modifiable parameters including filters, clusters, snapshots, and many more. Open this file and take a look. If this is your first time using PyMGal, your main focus should be on the "sim_dir" and "cat_dir" parameters. Make sure directories properly reflect the absolute path to your database. You can play around with other parameters, but the default values should be enough to get you started. The file can be found at your/path/to/pymgal/pymgal (i.e. the inner pymgal directory). 



#### Database configuration


- For a given cluster at a given redshift, you'll need the corresponding snap_XYZ simulation snapshot file. Once this is set up, the only thing missing is to find the coordinates which on which to center our projections. By default, PyMGal assumes you don't know these coordinates and therefore infers them from the halo catalogue files. To do this, make sure "custom_coords" is set to False in the config file and include the following files in your database. In addition to the simulation snapshot files, you'll need (2) one .AHF_halos file for each object and snapshot, (3) one .AHF_mtree_idx file for each object and snapshot, and (5) a single center-cluster.txt file which contains the halo IDs for each object at redshift 0.

<br />
```bash    
         ├── sim_dir
            ├── NewMDCLUSTER_0001
                ├── snap_128.hdf5
                ├── snap_127.hdf5
                ├── more snapshot files 
            ├── more NewMDCLUSTER_XXXX folders containing snapshot files
        
        ├── cat_dir 
            ├── GIZMO_R200c_snaps_128-center-cluster.txt
            ├── NewMDCLUSTER_0001
                ├── GIZMO-NewMDCLUSTER_0001.snap_128.z0.000.AHF_halos
                ├── GIZMO-NewMDCLUSTER_0001.snap_128.z0.000.AHF_mtree_idx
                ├── AHF_halos, and AHF_mtree_idx files for other snap numbers...
            ├── more NewMDCLUSTER_XXXX folders containing .AHF_halos and .AHF_mtree_idx files
               
```
<br />

- If you'd rather use your own custom coordinates, set "custom_coords" to True in your config file. In this case, all you need in your catalogue directory is a file called "coords.txt". This file must be in json format and look like {"object_0001.snap_128":[x1, y1, z1, r1], "object_0001.snap_127":[x2, y2, z2, r2]}. This is a strict requirement, so make sure you obey the format.

<br />
```bash  
        ├── cat_dir 
            ├── coords.txt
```



## How do I run the code?

* Once everything is ready, you can begin generating mock observations. To do so, run the following at command line.

    *  **`python /your/path/to/pymgal/pymgal generate_mocks.py <output directory>`**


* You can pass optional command line arguments to temporarily overrule the config file. You can even pass your own config file path as an argument and use that one instead.  To get more information about command line arguments, you can run the following. <br>

    * **`python generate_mocks.py --help`**
<br />

## Who do I talk to?

*   Please report any issue to Weiguang Cui cuiweiguang@gmail.com.
*   Or report a bug through issues.

## Acknowledgement

*  This package borrowed a lot things from ezgal (<http://www.baryons.org/ezgal/>). Please make your acknowledgement to their work when you use this package.
