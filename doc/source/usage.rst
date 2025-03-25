.. _usage:


Usage
==========



The MockObservation object
-------------

In most cases, the only API needed for PyMGal is the MockObservation object. MockObservation objects require two mandatory parameters: the path to your snapshot data and the coordinates and radius of the region you want to consider. If you don't know the coordinates of your object, you'll probably need to obtain some catalogue data.

When defining the path to your snapshot data, there are two options. If your simulation data is contained in a single file (e.g. snap_XYZ.hdf5), you can simply include the path to this file. If it is split between several snapshot files (e.g. snap_XYZ.0.hdf5, snap_XYZ.1.hdf5), you should specify the directory where these files are contained. If you do this, make sure to include every piece of the snapshot in the directory and no unrelated files. 

Once you initialize the object, you can generate projection files using the project() function. You may also optionally get the magnitudes of stellar particles with the get_mags() function, spectral energy distributions (SEDs) with the get_seds() function, or the physical properties (e.g. positions, masses, ages, metallcities) with the get_simdata() function. 
Here is a sample to get you started. 

.. code-block:: python

   from pymgal import MockObservation

   obs = MockObservation("/path/to/snapshot", [x_c, y_c, z_c, r])   
   obs.params["out_val"] = "luminosity"         # This is how you modify parameters
   obs.get_mags()                               # If you want to get the magnitudes of stellar particles in different filters
   obs.get_seds()                               # If you want to get the spectral energy distributions (SEDs) of stellar particles
   obs.get_simdata()                            # If you want to get the positions, masses, ages, and metallicities of stellar particles
   obs.project("/path/to/output")               # To generate and save the mock observation  



Note that after the SEDs have been calculated, you can modify many of these parameters (e.g. obs.params["proj_vecs"] = "y") without needing to recalculate the SEDs. This can save compute time. Note that this only works for parameters that are not directly involved in calculating SEDs such as projection effects.

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
               "rest_frame": True,
               "AR": 1.2,
               "npx": 512,
               "z_obs": 0.1,
               "ksmooth": 100,
               "g_soft": None,
               "p_thick": None,
               "add_spec": False,
               "spec_res": None,
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
    **Type: pymgal.SSP_model object.** The stellar population model you want to assume. PyMGal supports various types of models depending on your scientific goals. For more details, see :ref:`SSP Models <ssp_models>`.

- **dustf**:  
    **Type: pymgal.dusts object.** dustsThe dust attenuation function used in the model. Options include no dust, charlot_fall (Charlot and Fall 2000), or calzetti (Calzetti et al. 2000). You can also define a custom function within the dusts.py file if needed. 

- **filters**:  
    **Type: list.** The telescope filters you want to mimic in your mock observations. For more details, see :ref:`Filters <filters>`

- **out_val**:  
    **Type: string or None.** The units for the output data. Options include `"luminosity`" (erg/s), `"Lsun`" (solar luminosities), `"flux"` (erg/s/cm^2), `"jy"` (Jansky), `"Fv"` (erg/s/cm^2/Hz), `"Fl"` (erg/s/cm^2/angstrom), or `"magnitude"`. This is case-insensitive.

- **mag_type**:  
    **Type: string.** If `out_val` is set to `"magnitude"`, this parameter specifies the magnitude type. Options are `"AB"`, `"vega"`, `"solar"`, `"AB_app"`, `"vega_app"`, or `"solar_app"`. If `out_val` is not `"magnitude"`, this parameter has no effect.
    
- **proj_vecs**:  
    **Type: list.** A list of projection vectors. You can specify principal axes (i.e. "x", "y", or "z") or provide custom vectors in Cartesian coordinates [x, y, z]. Example usage: ["x", "y", "z", [0, -1, 0]]

- **rest_frame**:  
    **Type: bool.** If set to `True`, the magnitudes will be computed in the rest frame. Otherwise, they will be in the observer's frame.

- **AR**:  
    **Type: float.** The angular resolution in arcseconds. If set to None, it is automatically calculated. If both `AR` and `npx` are `"auto"`, `npx` defaults to 512.
    
- **npx**:  
    **Type: int.** The number of pixels in the output image. You can also set this to `"auto"`, which will automatically decide the pixel number to include all particles.
    
- **z_obs**:  
    **Type: float or None.** The redshift of the observation from the observer's point of view. This parameter affects only the apparent distance, not age or evolution. If set to None, it defaults to `max(0.05, sim z)`.

- **ksmooth**:  
    **Type: int.**  Must be non-negative. A smoothing parameter used in kNN Gaussian smoothing. The larger the `ksmooth` value, the smoother the results. Setting ksmooth=0 will turn the smoothing feature off. 

- **g_soft**:  
    **Type: int or None.** The gravitational softening length of the simulation in comoving kpc/h. This limits smoothing for mass/age/metal maps. If set to None, mass/age/metal are smoothed similarly to light.

- **p_thick**:  
    **Type: int.** The thickness cut (in kpc/h) along the projection direction. This cut is applied as `[center-p_thick, center+p_thick]`. If set to `null`, no cut is applied, and all data is used.
    
- **add_spec**:  
    **Type: bool.** Do you want to output the spectrum of your observation for your choice of axes? If set to True, the spectrum will be written to a fits file. 
    
- **spec_res**:  
    **Type float.** If you chose `add_spec = True`, you can modify the spectral resolution here. Accepted values: float: (0, 1], or an array in wavelength for sampling with unit of Hertz. If `add_spec = False`, this does nothing.

- **ncpu**:  
    **Type: int.** The number of CPUs used in parallel processing.

- **noise**:  
    **Type: list or None.** The noise level of Gaussian noise for your observations as a limiting AB mag at some specified signal to noise ratio within some specified circular aperture. Format: [mag_lim, SNR, r_aperture].

- **outmas**:  
    **Type: bool.** If `True`, the mass map corresponding to your data will be output.

- **outage**:  
    **Type: bool.** If `True`, the age map corresponding to your data will be output.

- **outmet**:  
    **Type: bool.** If `True`, the metallicity map corresponding to your data will be output.
    
- **quiet**:  
    **Type: bool.** If `True`, the print statements displaying the code's progress will be silenced
    
