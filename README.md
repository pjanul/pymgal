# README

PyMGal is a package that uses simple stellar synthesis models to generate mock observations of galaxies from hydrodynamical simulations. If you want a more detailed explanation of PyMGal than can be provided in a short readme file, see the documentation at [https://pymgal.readthedocs.io](https://pymgal.readthedocs.io).


If you're viewing this from the Github page, you may wish to view PyMGal's official page at the [Python Package Index (PyPI)](https://pypi.org/project/pymgal/).

If you're viewing this from the PyPI page, you may want to see the full [Github repository](https://github.com/pjanul/pymgal).


![pymgal_demo](https://github.com/user-attachments/assets/4e1a7977-c389-41a6-a644-edadb00f03a7)


Installation 
============

Since PyMGal is registered with PyPI, you can simply pip install the latest stable version. Make sure you have all the necessary prerequisites from the requirements.txt file beforehand. 

```python
pip install pymgal
```
 
Usage
============

In most cases, the only API needed for PyMGal is the MockObservation object. MockObservation objects require two mandatory parameters: the path to your snapshot file and the coordinates and radius of the region you want to consider. If you don't know the coordinates of your object, you'll probably need to obtain some catalogue data.

Once you initialize the object, you can generate projection files using the project() function. You may also optionally get the magnitudes of stellar particles with the get_mags() function, spectral energy distributions (SEDs) with the get_seds() function, or the physical properties (e.g. positions, masses, ages, metallcities) with the get_simdata() function. 
Here is a sample to get you started. 

```python
from pymgal import MockObservation

obs = MockObservation("/path/to/snapshot", [x_c, y_c, z_c, r])   
obs.params["out_val"] = "luminosity"         # This is how you modify parameters
obs.get_mags()                               # If you want to get the magnitudes of stellar particles in different filters
obs.get_seds()                               # If you want to get the spectral energy distributions (SEDs) of stellar particles
obs.get_simdata()                            # If you want to get the positions, masses, ages, and metallicities of stellar particles
obs.project("/path/to/output")               # To generate and save the mock observation  
```

If all goes well, you should see at least one newly formed snap_{XYZ}-{proj_angle}-{filter}.fits file in your output directory. 


Modifiable parameters
-------------

There are many different parameters you can modify for your magnitude calculations and your projections. Here is a list of them. For more information, see the [documentation website](https://pymgal.readthedocs.io). 

```python
def __init__(self, sim_file, coords, params=None):
        # Default parameter values
        defaults = {
            "model": SSP_models('bc03', IMF='chab', has_masses=True),
            "dustf": None,
            "filters": ["sdss_r"],
            "out_val": "flux",
            "mag_type": "AB",
            "proj_vecs": "z",
            "rest_frame": True,
            "AR": 1.2,
            "npx": 512,
            "z_obs": None,
            "ksmooth": 100,
            "g_soft": None,
            "psf": None,
            "add_spec": False,
            "spec_res": None,
            "p_thick": None,
            "ncpu": 16,
            "noise": None,
            "outmas": True,
            "outage": False,
            "outmet": False,
            "quiet": False
        }
```

Who do I talk to?
-----------

*   While PyMGal has been thoroughly tested by the developers, it is not impossible that bugs have find their way into the code.
*   If you find a bug, please report any issue to Weiguang Cui (cuiweiguang@gmail.com) or Patrick Janulewicz (patrick.janulewicz@gmail.com) or report it through issues.

Acknowledgement
----------
*  This package borrowed a lot things from EzGal (<http://www.baryons.org/ezgal/>). Please make your acknowledgement to their work when you use this package.

