# README


## What is this repository for?

* PyMGal (pronounced py-em-gal) is a package that uses simple stellar synthesis models to generate observed galaxies from hydrodynamical simulations.

## How do I get set up?

* Begin by downloading/cloning PyMGal onto your own device and placing it in your directory of choice. After downloading, you will need to install the necessary dependencies before you can run the code.

#### Dependencies 

* To install the necessary dependencies, simply enter your/path/to/pymgal (i.e. the outer PyMGal directory) and run the following at the command line.

    * `pip install -r requirements.txt`


#### The config.yaml file

* The config.yaml file contains modifiable parameters including the coordinates of your projection region, filters, SSP models, and many more. Open it and take a look. You can play around with these parameters, but as long as you have valid coordinates, the other default values should be enough to get you started. The file can be found at your/path/to/pymgal/pymgal/config.yaml (i.e. the inner PyMGal directory). 

* The full list of filters and SSP models can be found in the pymgal/filters and pymgal/models directories, respectively. 


## How do I run the code?

* Once everything is set up, you can begin generating mock observations. To do so, enter the inner PyMGal directory /your/path/to/pymgal/pymgal and run the following at command line.

    *  **`python generate_mocks.py <snapshot_file> <config_file> <output_dir>`**



* You can also pass optional command line arguments to temporarily overrule the config file. To get more information about command line arguments, you can run the following. <br>

    * **`python generate_mocks.py --help`**
    

#### Input:

* The path to your snapshot simulation file. Can be any flavour or Gadget or GIZMO. Can be formatted snap_XYZ or snap_XYZ.hdf5. 
* The path to your config.yaml file. You can use the built-in config.yaml file, or you can specify a different path and use that one. 
* The directory where you'd like to output your files.

#### Output:
* One FITS file for each selected projection angle and filter. File names will be formatted snap_{XYZ}-{proj_angle}-{filter}.fits. 

## What if I don't know the coordinates for my projections?

* In this case, you'll probably need halo catalogue data. Halo catalogues come in many formats including AHF (Amiga Halo Finder), FoF (Friends of Friends), Rockstar, and more. These catalogues will contain information regarding the physical positions and merger history of the particles in your simulation. You'll need to use these catalogues to obtain the physical coordinates of whatever object you'd like to project.

* ** Note: If you're working with data from The Three Hundred Project, we've included a script called pymgal/doc/the300_helper.py which helps you get positions from AHF halos. Open it and read the comments at the top for instructions. **

## Who do I talk to?

*   Please report any issue to Weiguang Cui (cuiweiguang@gmail.com) or Patrick Janulewicz (patrick.janulewicz@gmail.com).
*   Or report a bug through issues.

## Acknowledgement

*  This package borrowed a lot things from ezgal (<http://www.baryons.org/ezgal/>). Please make your acknowledgement to their work when you use this package.

