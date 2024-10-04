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

PyMGal supports various model types from different works. Below is a list of models that were created as a result of the EzGal package. They are BC03 from  Bruzual & Chalot (2003), M05 from Maraston (2005), CB07 from Charlot & Bruzual (2007), BaSTI from Percival et al. (2009), C09 from Conroy et al. (2009) and P2 for the PEGASE2 set from Fioc & Rocca-Volmerange (1997). 


For more details on these, consult the EzGal paper. If you'd like to download EzGal models, you can access them here: http://www.baryons.org/ezgal/download.php. Make sure to cite the EzGal authors and model authors.

.. list-table::
   :widths: 10 15 10
   :header-rows: 1

   * - Category
     - Available IMFs
     - Available CSPs
   * - P09
     - Kroupa
     - 
   * - BC03
     - Chabrier, Salpeter
     - Burst, constant, exponential
   * - Binary ised BC03 
     - Chabrier, Salpeter
     - 
   * - C09
     - Chabrier, Kroupa, Salpeter
     - Exponential
   * - CB07
     - Chabrier, Salpeter
     - Burst, constant, exponential
   * - M05
     - Kroupa, Salpeter
     - Exponential
   * - P2
     - Salpeter
     - 


 
  
.. _dust_funcs:

Dust attentuation functions
----------

By default, PyMGal doesn't account for dust attenuation. If you want to add the effect of dust, you can use either the dust function described in Charlot and Fall (2000) or Calzetti et al. (2000). 

If you want to code your own dust function, you should be able to add it to the dusts.py file and then call it when creating your MockObservation object.
 
 
   
