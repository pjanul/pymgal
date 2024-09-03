Parameters
==========

This document describes the various parameters used in PyMGal for generating optical mock observations. Each parameter plays a specific role in defining the characteristics of the simulation, projection, and output.

Coordinate Parameters
----------------------

- **coords**:  
    The coordinates `[x, y, z, r]` of the projection region. Here, `(x, y, z)` represents the center of the projection region, and `r` is its radius. It's recommended to overwrite this with the correct coordinates from a halo file.

Model and Filter Parameters
---------------------------

- **SSP_model**:  
    The Simple Stellar Population (SSP) model you wish to use. The full list of available models can be found in `pymgal//models`.

- **IMF**:  
    The Initial Mass Function (IMF) assumed by the SSP model. Available options are `"chab"` (Chabrier), `"krou"` (Kroupa), or `"salp"` (Salpeter).

- **dust_func**:  
    The dust attenuation function used in the model. Options include `"null"` (no dust), `"charlot_fall"` (Charlot and Fall 2000), or `"calzetti"` (Calzetti et al. 2000).

- **filters**:  
    The telescope filters you want to mimic in your mock observations. You can list multiple filters here.

- **out_val**:  
    The units for the output data. Options include `"flux"` (erg/s/cm^2), `"jy"` (Jansky), `"Fv"` (erg/s/cm^2/Hz), `"Fl"` (erg/s/cm^2/Ångstrom), or `"magnitude"`. This is case-insensitive.

- **mag_type**:  
    If `out_val` is set to `"magnitude"`, this parameter specifies the magnitude type. Options are `"AB"`, `"vega"`, `"solar"`, `"AB_app"`, `"vega_app"`, or `"solar_app"`. If `out_val` is not `"magnitude"`, this parameter has no effect.

Projection Parameters
----------------------

- **num_random_proj**:  
    The number of random projections you want to generate. Setting this to `0`, along with `proj_vecs = null` and `proj_angles = null`, will cause an error.

- **proj_vecs**:  
    A list of projection vectors. You can specify principal axes (`"x"`, `"y"`, `"z"`) or provide custom vectors in Cartesian coordinates `[x, y, z]`. This can also be set to `null`.

- **proj_angles**:  
    A list of angles `[alpha, beta, gamma]` (in degrees) used to rotate around the x, y, and z axes, respectively. This serves the same purpose as `proj_vecs`, so you can use either or both. This can also be set to `null`.

Miscellaneous Parameters
------------------------

- **rest_frame**:  
    If set to `True`, the results will be in the rest frame. Otherwise, they will be in the observer's frame.

- **AR**:  
    The angular resolution in arcseconds. If set to `null`, it is automatically calculated. If both `AR` and `npx` are `"auto"`, `npx` defaults to 512.

- **ksmooth**:  
    A smoothing parameter used in kNN Gaussian smoothing. The larger the `ksmooth` value, the smoother the results.

- **g_soft**:  
    The gravitational softening length of the simulation in comoving kpc/h. This limits smoothing for mass/age/metal maps. If set to `null`, mass/age/metal are smoothed similarly to light.

- **z_obs**:  
    The redshift of the observation from the observer's point of view. This parameter affects only the apparent distance, not age or evolution. If set to `null`, it defaults to `max(0.05, sim z)`.

- **npx**:  
    The number of pixels in the output image. You can also set this to `"auto"`, which will automatically decide the pixel number to include all particles.

- **zthick**:  
    The thickness cut (in kpc/h) along the projection direction. This cut is applied as `[center-zthick, center+zthick]`. If set to `null`, no cut is applied, and all data is used.

- **ncpu**:  
    The number of CPUs used in parallel processing.

- **outmas**:  
    If `True`, the mass map corresponding to your data will be output.

- **outage**:  
    If `True`, the age map corresponding to your data will be output.

- **outmet**:  
    If `True`, the metallicity map corresponding to your data will be output.
