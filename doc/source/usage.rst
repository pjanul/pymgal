.. _usage:


Usage
==========



The MockObservation object
-------------

In most cases, the only API needed for PyMGal is the MockObservation object. MockObservation objects require two mandatory parameters: the path to your snapshot file and the coordinates and radius of the region you want to consider. If you don't know the coordinates of your object, you'll probably need to obtain some catalogue data.

Once you initialize the object, you can calculate magnitudes of particles in your preferred output unit using the get_mags() function. You can also save projection files using the project() function. If you call project() before calling get_mags(), the magnitudes will automatically be calculated.

Here is a sample to get you started. If all goes well, you should see at least one newly formed snap_{XYZ}-{proj_angle}-{filter}.fits file in your output directory.

.. code-block:: python

   from pymgal import MockObservation

   obs = MockObservation("/path/to/snapshot", [x_c, y_c, z_c, r])   
   obs.params["out_val"] = "luminosity"
   obs.get_mags()
   obs.project("/path/to/output")


List of parameters
-------------

There are many different parameters you can modify for your magnitude calculations and your projections. You can find the full list below. 


.. code-block:: python

   class MockObservation(object):
       def __init__(self, sim_file, coords, args=None):
           # Default parameter values
           defaults = {
               "model": SSP_models('bc03', IMF='chab', has_masses=True),
               "dustf": None,
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
               "p_thick": None,
               "ncpu": 16,
               "noise": None,
               "outmas": True,
               "outage": False,
               "outmet": False
               "quiet": False
           }

Description of parameters
-----------

This document describes the various parameters used in PyMGal for generating optical mock observations. Each parameter plays a specific role in defining the characteristics of the simulation, projection, and output.



- **model**:  
    The stellar population model you want to assume. PyMGal supports various types of models. For more details, see :ref:`SSP Models <ssp_models>`.

- **dustf**:  
    The dust attenuation function used in the model. Options include no dust, "charlot_fall" (Charlot and Fall 2000), or "calzetti" (Calzetti et al. 2000). You can also define a custom function within the dusts.py file if needed.

- **filters**:  
    The telescope filters you want to mimic in your mock observations. For more details, see :ref:`Filters <filters>`

- **out_val**:  
    The units for the output data. Options include `"luminosity`" (erg/s), `"Lsun`" (solar luminosities), `"flux"` (erg/s/cm^2), `"jy"` (Jansky), `"Fv"` (erg/s/cm^2/Hz), `"Fl"` (erg/s/cm^2/angstrom), or `"magnitude"`. This is case-insensitive.

- **mag_type**:  
    If `out_val` is set to `"magnitude"`, this parameter specifies the magnitude type. Options are `"AB"`, `"vega"`, `"solar"`, `"AB_app"`, `"vega_app"`, or `"solar_app"`. If `out_val` is not `"magnitude"`, this parameter has no effect.
    
- **proj_vecs**:  
    A list of projection vectors. You can specify principal axes (i.e. "x", "y", or "z") or provide custom vectors in Cartesian coordinates [x, y, z].

- **proj_angs**:  
    A list of angles `[alpha, beta, gamma]` (in degrees) used to rotate around the x, y, and z axes, respectively. This serves the same purpose as `proj_vecs`, so you can use either or both.
    
- **proj_rand**:  
    The number of random projections you want to generate. Setting this to `0`, along with `proj_vecs = null` and `proj_angles = null`, will cause an error.

- **rest_frame**:  
    If set to `True`, the magnitudes will be computed in the rest frame. Otherwise, they will be in the observer's frame.

- **AR**:  
    The angular resolution in arcseconds. If set to `null`, it is automatically calculated. If both `AR` and `npx` are `"auto"`, `npx` defaults to 512.
    
- **npx**:  
    The number of pixels in the output image. You can also set this to `"auto"`, which will automatically decide the pixel number to include all particles.
    
- **z_obs**:  
    The redshift of the observation from the observer's point of view. This parameter affects only the apparent distance, not age or evolution. If set to `null`, it defaults to `max(0.05, sim z)`.

- **ksmooth**:  
    A smoothing parameter used in kNN Gaussian smoothing. The larger the `ksmooth` value, the smoother the results.

- **g_soft**:  
    The gravitational softening length of the simulation in comoving kpc/h. This limits smoothing for mass/age/metal maps. If set to `null`, mass/age/metal are smoothed similarly to light.

- **p_thick**:  
    The thickness cut (in kpc/h) along the projection direction. This cut is applied as `[center-p_thick, center+p_thick]`. If set to `null`, no cut is applied, and all data is used.

- **ncpu**:  
    The number of CPUs used in parallel processing.

- **noise**:  
    The noise level of Gaussian noise for your observations in AB mag / arcsec^2, which will be converted to your choice of out_val

- **outmas**:  
    If `True`, the mass map corresponding to your data will be output.

- **outage**:  
    If `True`, the age map corresponding to your data will be output.

- **outmet**:  
    If `True`, the metallicity map corresponding to your data will be output.
    
- **quiet**:  
    If `True`, the print statements displaying the code's progress will be silenced
    
