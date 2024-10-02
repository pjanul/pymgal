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
  
 
Usage
-------------

In most cases, the only API needed for PyMGal is the `MockObservation` object. Here's an example of how to import and generate the MockObservation object, as well as how to modify it's parameters. You can either modify the parameters as you create the object, or you can switch them after initializing. MockObservation requires two mandatory parameters: the path to your snapshot file and the coordimates + radius of the region you want to consider. Once you initialize the object, you can calculate magnitudes of particles in your preferred output unit using the get_mags() function. You can also save projection files using the project() function, provided you've given the function your output directoy. Here is a sample to get you started.

.. code-block:: python

   from pymgal import MockObservation

   obs = MockObservation("/path/to/snapshot", [x_c, y_c, z_c, r])   
   obs.params["out_val"] = "luminosity"
   obs.get_mags()
   obs.project("/path/to/output")


Modifiable parameters
-------------

There are many different parameters you can modify for your magnitude calculations and your projections. Here is a list of them. For more information, see the :ref:Parameters <parameters> page.


.. code-block:: python

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