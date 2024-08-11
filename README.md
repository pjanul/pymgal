# README


## What is this repository for?

* PyMGal (pronounced py-em-gal) is a package that uses simple stellar synthesis models to generate observed galaxies from hydrodynamical simulations.

## How do I get set up?

* Begin by downloading/cloning PyMGal onto your own device and placing it in your directory of choice. After downloading, you will need to install the necessary dependencies and configure your database before you can run the code.

#### Dependencies 

* To install the necessary dependencies, simply enter the outer pymgal directory and run the following at the command line.

    * `pip install -r requirements.txt`


#### The config.yaml file

* The config.yaml file contains modifiable parameters including filters, region IDs, snapshots, and many more. Open it and take a look. To get started, your main focus should be on the **`sim_dir`** and  **`coords_filepath`** parameters. We will provide instructions for setting these up. You can play around with other parameters, but the default values should be enough to get you started. The file can be found at your/path/to/pymgal/pymgal/config.yaml (i.e. the inner pymgal directory). 



#### Database configuration

* For a given region at a given redshift, you'll need the corresponding simulation snapshot file, which will either look like snap_XXX or snap_XXX.hdf5, depending on the file format. 

<br />
```bash    
         ├── sim_dir
            ├── RegionPrefix_0001
                ├── snap_XXX.hdf5
                ├── snap_XXY.hdf5
                ├── more snapshot files 
            ├── RegionPrefix_0002
                ├── more snapshot files     
```
<br />

* Now all you need is a file containing the coordinates which represent the centre and radius of your projection region. This file must be in json format and look like `{"RegionPrefix_0001.snap_XXX":[x1, y1, z1, r1], "RegionPrefix_0001.snap_XXY":[x2, y2, z2, r2], etc.}`. The formatting is important for the code to run properly, so make sure to follow this convention.

* Make sure to set **`sim_dir: path/to/sim_dir`** and  **`coords_filepath: path/to/coords.txt`** are updated to reflect these directories. 


## How do I run the code?

* Once everything is ready, you can begin generating mock observations. To do so, run the following at command line.

    *  **`python /your/path/to/pymgal/pymgal generate_mocks.py <output directory>`**


* You can pass optional command line arguments to temporarily overrule the config file. You can even pass your own config file path as an argument and use that one instead.  To get more information about command line arguments, you can run the following. <br>

    * **`python generate_mocks.py --help`**
<br />

## What if I don't know the coordinates for my projections?

* In this case, you'll probably need halo catalogue data. Halo catalogues come in many formats including AHF (Amiga Halo Finder), Friends of Friends (FoF), Rockstar, and more. These catalogues will contain information regarding the physical positions and merger history of the particles in your simulation. You'll need to use these catalogues to obtain the physical coordinates of whatever object you'd like to project.

    * **`Note: If you're working with data from The Three Hundred Project, contact Patrick to get the relevant coords.txt file and/or the script to generate it. **

## Who do I talk to?

*   Please report any issue to Weiguang Cui (cuiweiguang@gmail.com) or Patrick Janulewicz (patrick.janulewicz@mail.mcgill.ca).
*   Or report a bug through issues.

## Acknowledgement

*  This package borrowed a lot things from ezgal (<http://www.baryons.org/ezgal/>). Please make your acknowledgement to their work when you use this package.
