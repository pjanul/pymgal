
Models and outputs
==========

.. _ssp_models:
SSP Models
----------

In PyMGal, you can choose between different Simple Stellar Population (SSP) models and Initial Mass Functions (IMFs). These models are used to simulate the spectral energy distribution (SED) of stellar populations in galaxies. You can choose the a model based on your scientific goals, or simply use the default model.

Available SSP Models:
 * BC03 from  Bruzual & Chalot (2003)
 * M05 from Maraston (2005)
 * CB07 from Charlot & Bruzual (2007)
 * BaSTI from Percival et al. (2009)
 * C09 from Conroy et al. (2009)
 
.. _imf:
Initial Mass Functions
----------

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