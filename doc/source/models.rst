.. _ssp_models:
Models and outputs
==========

.. _ssp_models_object:

The SSP_models object
----------

If you want to change the stellar population model used by the program, you'll probably want to call the SSP_models object. 

You can use any of PyMGal's available models, or you can specify your own by adding a .txt file and specifying is_ascii=True in the object initialization.

.. code-block:: python

   from pymgal import SSP_models, MockObservation
   
   model_type = "bc03"
   model = SSP_models(model_type, IMF="chab")
   obs = MockObservation("/path/to/snapshot", [x_c, y_c, z_c, r], params = {"model": model})
   
.. _ssp_models_params:

SSP_model parameters
----------

The code below shows the full list of options for the SSP_models class, as well as how you can print its docstring for more details. 


.. code-block:: python

   model = SSP_models(model_file, IMF="chab", metal=[list], is_ised=False, is_fits=False,
   is_ascii=False, has_masses=False, units='a', age_units='gyrs', nsample=None, quiet=True, model_dir=None)
   
   print(SSP_models.__doc__)

The most important are the model_file and the IMF. If you want to only look at some metallicities, specify them with the metal parameter. If you want to specify your own model in a .txt file, make sure to set is_ascii=True. You may also want to specify your own custom model directory by setting model_dir="/path/to/your/models"


.. _avail_models:
Available models
----------
Here These models are used to simulate the spectral energy distribution (SED) of stellar populations in galaxies. You can choose the a model based on your scientific goals, or simply use the default model.

Available SSP Models:
 * BC03 from  Bruzual & Chalot (2003)
 * M05 from Maraston (2005)
 * CB07 from Charlot & Bruzual (2007)
 * BaSTI from Percival et al. (2009)
 * C09 from Conroy et al. (2009)
 
.. _imf:
Initial Mass Functions
----------

Note that not all models support all three IMFs.

Available IMFs:
 * Chabrier 
 * Kroupa
 * Salpeter
 
  
.. _dust_funcs:
Dust attentuation functions
----------

Available dust functions:
 * None
 * Charlot and Fall (2000)
 * Calzetti et al. (2000)
 
 
   
.. _out_vals:
Output units
----------

Output units:
 * Luminosity (erg/s)
 * Flux (erg/s/cm^2)
 * Spectral flux density
    * F_nu (erg/s/cm^2/Hz)
    * F_lambda (erg/s/cm^2/Angstrom)
    * Jansky 
 * Magnitude (either absolute or apparent)
    * AB
    * Solar
    * Vega