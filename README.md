# README


What is this repository for?
=============

* PyMGal is a package that uses simple stellar synthesis models to generate observed galaxies from hydrodynamical simulations.

Documentation
=============

If you want a more detailed explanation of PyMGal than can be provided in a short readme file, see the documentation at [https://pymgal.readthedocs.io](https://pymgal.readthedocs.io).

Installation 
============

Installing developer version
-------------
To install the latest version, you can clone the repository with git. 

  * git clone https://bitbucket.org/pjanul/pymgal
  
Prerequisites
-------------

To install the necessary dependencies, simply enter your/path/to/pymgal (i.e. the outer PyMGal directory) and run the following at the command line.

  * pip install -r requirements.txt
  
 
Usage
============

In most cases, the only API needed for PyMGal is the MockObservation object. MockObservation objects require two mandatory parameters: the path to your snapshot file and the coordimates + radius of the region you want to consider. If you don't know the coordinates of your object, you'll probably need to obtain some catalogue data.

Once you initialize the object, you can calculate magnitudes of particles in your preferred output unit using the get_mags() function. You can also save projection files using the project() function. If you call project() before calling get_mags(), the magnitudes will automatically be calculated.

Here is a sample to get you started. If all goes well, you should see at least one newly formed snap_{XYZ}-{proj_angle}-{filter}.fits file in your output directory.



```python
from pymgal import MockObservation

obs = MockObservation("/path/to/snapshot", [x_c, y_c, z_c, r])   
obs.params["out_val"] = "luminosity"
obs.get_mags()
obs.project("/path/to/output")
```


Modifiable parameters
-------------

There are many different parameters you can modify for your magnitude calculations and your projections. Here is a list of them. For more information, see documentation website here: [https://pymgal.readthedocs.io](https://pymgal.readthedocs.io).

```python 
class MockObservation(object):
    def __init__(self, sim_file, coords, args=None):
        # Default parameter values
        defaults = {
            "model": "bc03",
            "imf": "chab",
            "dustf": None,
            "custom_model": None,
            "filters": ["sdss_r"],
            "out_val": "flux",
            "mag_type": "AB",
            "proj_vecs": "z",
            "proj_angs": None,
            "proj_rand": 0,
            "rest_frame": True,
            "AR": 1.2,
            "npx": 512,
            "z_obs": 0.1,
            "ksmooth": 100,
            "g_soft": None,
            "thickness": None,
            "ncpu": 16,
            "noise": None,
            "outmas": True,
            "outage": False,
            "outmet": False
        }
```

What if I don't know the coordinates for my projections?
----------

* In this case, you'll probably need halo catalogue data. Halo catalogues come in many formats including AHF (Amiga Halo Finder), FoF (Friends of Friends), Rockstar, and more. These catalogues will contain information regarding the physical positions and merger history of the particles in your simulation. You'll need to use these catalogues to obtain the physical coordinates of whatever object you'd like to project.

* ** Note: If you're working with data from The Three Hundred Project, we've included a script called pymgal/doc/the300_helper.py which helps you get positions from AHF halos. Open it and read the comments at the top for instructions. **

Who do I talk to?
-----------

*   Please report any issue to Weiguang Cui (cuiweiguang@gmail.com) or Patrick Janulewicz (patrick.janulewicz@gmail.com).
*   Or report a bug through issues.

Acknowledgement
----------
*  This package borrowed a lot things from ezgal (<http://www.baryons.org/ezgal/>). Please make your acknowledgement to their work when you use this package.

