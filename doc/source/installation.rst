Installation
==================

This guide will walk you through the necessary prerequisites, installation steps, and how to run the code for PyMGal.

Installing stable version
-------------
We are working on registering PyMGal with the Python Package Index (PyPI). Once this is done, PyMGal will be installable with pip. Until then, please install the developer version.

Installing developer version
-------------
To install the latest version, you can clone the repository with git. 

  * git clone https://bitbucket.org/pjanul/pymgal
  
Prerequisites
-------------

To install the necessary dependencies, simply enter your/path/to/pymgal (i.e. the outer PyMGal directory) and run the following at the command line.

  * pip install -r requirements.txt
  
  
The config.yaml file
-------------
The config.yaml file contains modifiable parameters including the coordinates of your projection region, filters, SSP models, and many more. Open it and take a look. You can play around with these parameters, but as long as you have valid coordinates, the other default values should be enough to get you started. The file can be found at your/path/to/pymgal/pymgal/config.yaml (i.e. the inner PyMGal directory). For more details on this, see the Parameters page.


Running the code
-------------

* Once everything is set up, you can begin generating mock observations. To do so, enter the inner PyMGal directory /your/path/to/pymgal/pymgal and run the following at command line.

    *  **`python generate_mocks.py <snapshot_file> <config_file> <output_dir>`**



* You can also pass optional command line arguments to temporarily overrule the config file. To get more information about command line arguments, you can run the following. <br>

    * **`python generate_mocks.py --help`**
    

* Input:

   * The path to your snapshot simulation file. Can be any flavour or Gadget or GIZMO. Can be formatted snap_XYZ or snap_XYZ.hdf5. 
   * The path to your config.yaml file. You can use the built-in config.yaml file, or you can specify a different path and use that one. 
   * The directory where you'd like to output your files.

* Output:
   * One FITS file for each selected projection angle and filter. File names will be formatted snap_{XYZ}-{proj_angle}-{filter}.fits. 
