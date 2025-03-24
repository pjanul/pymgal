import os
import numpy as np
import astropy.io.fits as pyfits
from pymgal import utils
from scipy.interpolate import interp1d
from scipy.integrate import simpson
# from astropy.cosmology import FlatLambdaCDM



class filters(object):
    r""" load filters for telescopes.
    filter = SSP_model(path, name)

    Parameters
    ----------
    path : The folder where these filters are saved. Default: ''.
           The program will take the buildin filter folder or the specified
           filter environment.
    name : The name/s of the filter. Default: ''. The program will not load
           filters. But you can add it later by filter.add_filter(). But the
           filter folder can not be changed in that function.

    If the file path is not found then program will search for it in the directory
    specified by the ``data/filters`` and ``ezsps_FILTERS`` environment variable,
    at last the directory in the program module directory /filters/.
    """

    def __init__(self, f_path='', f_name='', units='a'):
        # clear filter list etc
        self.f_vs = {}
        self.f_ls = {}
        self.f_tran = {}
        self.npts = {}
        self.filter_order = []
        self.nfilters = 0              # number of filters
        self.current_filter = -1       # counter for iterator
        self.ab_flux = {}
        self.vega_mag = {}
        self.solar_mag = {}
        self.noises = {}
        self.spectrum = {}
        
        # tolerance for determining whether a given zf matches a stored zf
        # the tolerance is typical set by ezsps after creating a new astro filter
        # but it is also defined here to have a default value
        self.tol = 1e-8
        self.ab_source_flux = 3.631e-20  # flux of a zero mag ab source

        # save path to data folder: module directory/data
        self.data_dir = os.path.dirname(os.path.realpath(__file__)) + '/'
        # make sure paths end with a slash
        if self.data_dir[-1] != os.sep:
            self.data_dir += os.sep

        # how about path to filter and model directories?
        self.filter_dir = False
        if not f_path:
            if 'ezsps_filters' in os.environ:
                self.filter_dir = os.environ['ezsps_filters']
            elif 'ezsps_FILTERS' in os.environ:
                self.filter_dir = os.environ['ezsps_FILTERS']
            else:
                self.filter_dir = '%sfilters/' % self.data_dir
        else:
            self.filter_dir = f_path

        if self.filter_dir and self.filter_dir[-1] != os.sep:
            self.filter_dir += os.sep

        # attempt to load the vega spectrum
        vega_file = '%srefs/vega.fits' % self.data_dir
        if os.path.isfile(vega_file):
            fits = pyfits.open(vega_file)
            self.vega = np.column_stack(
                (fits[1].data.field('freq'), fits[1].data.field('flux')))
            self.has_vega = True
        else:
            self.vega = np.array([])
            self.has_vega = False

        # attempt to load the solar spectrum
        solar_file = '%srefs/solar.fits' % self.data_dir
        if os.path.isfile(solar_file):
            fits = pyfits.open(solar_file)
            self.solar = np.column_stack(
                (fits[1].data.field('freq'), fits[1].data.field('flux')))
            self.has_solar = True
        else:
            self.solar = np.array([])
            self.has_solar = False

        if f_name:
            self.add_filter(f_name, units=units)

    # #####################
    # #  return iterator  #
    # #####################
    # def __iter__(self):
    #     self.current_filter = -1
    #     return self
    #
    # #########################
    # #  next() for iterator  #
    # #########################
    # def next(self):
    #
    #     self.current_filter += 1
    #     if self.current_filter == len(self.filter_order):
    #         raise StopIteration
    #
    #     filt = self.filter_order[self.current_filter]
    #     return (filt, self.filters[filt])

    def add_filter(self, name, units='a', grid=True):
        r"""
        ezsps.add_filter(name, units='a', grid=True)

        :param name: The name to store the filter as. String or list of strings
        :param units: The length units for the wavelengths in the file. String
        :param grid: Whether or not to calculate evolution information when first added

        Add a filter for calculating models.
        Specify the name of the file containing the filter transmission curve.

        The filter file should have two columns (wavelength,transmission).
        Wavelengths are expected to be in angstroms unless specified otherwise
        with ``units``.
        See :func:`ezsps.utils.to_meters` for list of available units.

        Specify a name to refer to the filter as later.
        If no name is specified, the filename is used (excluding path information
        and extension)
        If a filter already exists with that name, the previous filter will be replaced.

        If grid is True, then models will be generated for this filter at all
        set formation redshifts.

        You can pass a numpy array directly, instead of a file,
        but if you do this you need to specify the name.
        """

        if isinstance(name, type('')):  # single filter
            if not self.filter_order.count(name):
                self.filter_order.append(name)
                self.nfilters += 1
                self._load_filter(name, units=units)
        elif isinstance(name, type([])):  # multiple filters
            for i in name:
                if not self.filter_order.count(i):
                    self.filter_order.append(i)
                    self.nfilters += 1
                    self._load_filter(i, units=units)
        else:
            raise ValueError(
                'You need to pass a filter name or a list of filter names!')

        # self.filters[name].tol = self.tol
        # store its name in self.filter_order

        # if grid:
        #     self._grid_filters(name)
    def _load_filter(self, fname, units='a'):
        filters = utils.rascii(self.filter_dir + fname)

        # calculate wavelengths in both angstroms and hertz
        units = units.lower()
        if units == 'hz':
            self.f_vs[fname] = filters[:, 0]
            self.f_ls[fname] = utils.to_lambda(self.f_vs[fname], units='a')
            self.f_tran[fname] = filters[:, 1]
        else:
            self.f_vs[fname] = utils.to_hertz(filters[:, 0], units='a')[::-1]
            self.f_ls[fname] = utils.convert_length(filters[:, 0], incoming=units, outgoing='a')[::-1]
            self.f_tran[fname] = filters[:, 1][::-1] # need to convert from small to large in frequency the same as SSP model

        self.npts[fname] = filters[:, 0].size

        # normalization for calculating ab mags for this filter
        self.ab_flux[fname] = self.ab_source_flux * \
            simpson(self.f_tran[fname] / self.f_vs[fname], x=self.f_vs[fname])

        # store the cosmology object if passed
        # if cosmology is not None:
        #     self.cosmo = cosmology

        # calculate ab-to-vega conversion if vega spectrum was passed
        if self.has_vega:
            if self.f_vs[fname].min() < self.vega[:, 0].min() or \
               self.f_vs[fname].max() > self.vega[:, 0].max():
                raise ValueError(
                    'The filter frequency is out of Vega frequency!')
            self.vega_mag[fname] = -1.0 * self._VS_mag(fname, vega=True)
            # self.calc_mag(self.vega[:, 0],
            #               self.vega[:, 1], 0, fn=fname)[fname]
            # self.vega_flux = self.ab_source_flux / 10.0**(-0.4 * self.vega_mag)

        # calculate solar magnitude if solar spectrum was passed
        if self.has_solar:
            # does the solar spectrum extend out to this filter?
            if self.f_vs[fname].min() < self.solar[:, 0].min() or \
               self.f_vs[fname].max() > self.solar[:, 0].max():
                raise ValueError(
                    'The filter frequency is out of Solar frequency!')
            self.solar_mag[fname] = self._VS_mag(fname)
            # self.calc_mag(self.solar[:, 0], self.solar[:, 1], 0, fn=fname)[fname]

    ##############
    #  calc vega solar mag  #
    ##############
    def _VS_mag(self, fn, vega=False):
        r""" Default calculate solar mag"""

        if vega:
            vs = self.vega[:, 0]
            se = self.vega[:, 1]
        else:
            vs = self.solar[:, 0]
            se = self.solar[:, 1]
        interp = interp1d(vs, se, axis=0, bounds_error=False, fill_value="extrapolate")
        sed_flux = simpson(interp(self.f_vs[fn]).T * self.f_tran[fn] / self.f_vs[fn], x=self.f_vs[fn])
        return -2.5 * np.log10(sed_flux / self.ab_flux[fn])
    
    ##############
    #  calc energy in filter  #
    ##############
    def calc_energy(self, sspmod, simd, fn=None, dust_func=None, unit='flux',
                    apparent=False, vega=False, solar=False, rest_frame=False, noise=None, redshift=0.05, outspec_res=None):
        r"""
        mag = pymgal.calc_energy(self, sspmod, simd, fn=None, dust_func=None, unit='flux',
                                 apparent=False, vega=False, solar=False, rest_frame=False):

        sspmod      : loaded ssp models.
        simd        : loaded simulation data from load_data
        fn          : the name of filter/s. string, list of strings, or None
                        By default (None), all loaded filters will be included in the calculation
        dust_func   : The function for dust attenuetion.
        unit        : the returned value in the filter bands. Can be 'luminosity', 'flux', 'magnitude', 'Fv', 'Fl', or 'Jy'
        apparent    : If you need apparent magnitude, set True. Default False
        vega        : return AB magnitude in VEGA.
        solar       : return AB magnitude in solar.
        rest_frame  : Do you want the resutls in rest_frame? default: False, in observer's frame.
                        We always assume the observer is at z = 0. If true, we force everything in rest_frame.
        redshift    : The redshift specified by the user.
        outspec_res : The energy resolution for outputting spectrum. Default: None, the model's finest resolution is output
                        Accepted values: float: (0, 1], or an array in wavelength for sampling with unit of Hertz.
        returns     : The brightness in the units specified by the user. 
                      Magnitude by default is AB, can be apparent if apparent=True, in vega mags if vega=True or in solar mags in solar=True

        
        Available output units are (case insensitive):

        ========== ===============================
        Name       Units
        ========== ===============================
        Jy         Jansky
        Fv         ergs/s/cm^2/Hz
        Fl         ergs/s/cm^2/Angstrom
        Flux       ergs/s/cm^2
        Luminosity ergs/s
        Magnitude  AB solar, or vega in absolute or apparent
        ========== ================================

        Example:
            mag = pymgal.calc_mag(ssp, simd, z, fn=None)
        
        Notes:
        Calculate the flux, spectral flux density, luminosity or magnitude of the given sed.
        """

        # make sure an acceptable number of sed points actually go through the
        # filter...
        if isinstance(fn, type('')):  # single filter will be used
            fn = [fn]
        else:
            if fn is None:
                fn = self.filter_order
            else:
                if not isinstance(fn, type([])):
                    raise ValueError("Incorrected filter name ! ", fn)

        # if not redshift:
        #     redshift = simd.z
        # if redshift < 0:
        #     redshift = 0

        if apparent:
            app = 5. * np.log10(simd.cosmology.luminosity_distance(redshift).value
                                * 1.0e5) if redshift > 0 else -np.inf
        else:
            app = 0.0


        vs = sspmod.vs[sspmod.met_name[0]]
        mag = {}
        
        vsn, sedn = sspmod.get_seds(simd, rest_frame=rest_frame, dust_func=dust_func)
        
        units=unit.lower() # case insensitive
        
        # If we want to get rid of the Hz dependence, convert to erg/s
        if units in {'luminosity', 'flux', 'lsun'}:
            sedn *= sspmod.vs[sspmod.met_name[0]].reshape(sspmod.nvs[0], 1) 
        elif units == "jy":
            sedn *= 1e23
        elif units == "fl":
            sedn *= utils.convert_length(utils.c, outgoing='a') / sspmod.ls[sspmod.met_name[0]].reshape(sspmod.nls[0], 1)**2 

        vsn = np.asarray(vsn, dtype='<f8')
        sedn = np.asarray(sedn, dtype='<f8')
        interp = utils.numba_interp1d(vsn, sedn.T)
        #end_time = time.time()
        #print("filters_time: ", end_time - start_time)
        
        for i in fn:
            if vega:
                to_vega = self.vega_mag[i]
            else:
                to_vega = 0.0
            if solar:
                to_solar = self.solar_mag[i]
            else:
                to_solar = 0.0
 

            d_L = simd.cosmology.luminosity_distance(redshift).to('cm').value  # convert user-defined redshift to a luminosity distance in cm (for luminosity -> flux)
            L_sun = 3.826e33 # solar luminosity
        
            # Handle conversions for Gaussian noise
            lambda_c = simpson(self.f_ls[i] * self.f_tran[i], x=self.f_ls[i]) / simpson(self.f_tran[i], x=self.f_ls[i])
            nu_c = simpson(self.f_vs[i] * self.f_tran[i], x=self.f_vs[i]) / simpson(self.f_tran[i], x=self.f_vs[i])
            #noise_mag = None                                         # Define noise_mag so we don't call it later without initiating it first
            if noise is not None:
                noise_mag = noise[0] + to_vega + to_solar + app if units == "magnitude" else 10**(-0.4*(noise[0]+48.6))  #output magnitude if needed, otherwise convert to Fv. Recall noise=[mag_lim, SNR, r_aperture]
                noise_mag = (
                    noise_mag * nu_c * (4.0 * np.pi * d_L**2.0) if units == "luminosity" else
                    noise_mag * nu_c * (4.0 * np.pi * d_L**2.0) / L_sun if units == "lsun" else
                    noise_mag * nu_c if units == "flux" else
                    noise_mag * 10**23 if units == "jy" else
                    noise_mag * utils.convert_length(utils.c, outgoing='a') / lambda_c**2 if units == "fl"
                    else noise_mag         # if units are "fv"
                )
                self.noises[i] = [noise_mag, noise[1], noise[2]] # Create a dictionary where the keys are filter names and values are [mag_lim, SNR, r_aperture]
            else: self.noises = None
            
            
            c = ((vsn > self.f_vs[i].min()) & (vsn < self.f_vs[i].max())).sum()
            if c < 3:
                print("Warning, wavelength range from SSP models is too near filter,",
                      " magnitude is assigned nan")
                mag[i] = np.nan
            # and that the SED actually covers the whole filter
            if vsn.min() > self.f_vs[i].min() or vsn.max() < self.f_vs[i].max():
                print("Warning, wavelength range from SSP models is outside of filter,",
                      " magnitude is assigned nan")
                mag[i] = np.nan
            #print("nu_c", nu_c)
            # In units of solar luminosity? Verify this
            mag[i] = simpson(interp(self.f_vs[i]) * self.f_tran[i] / self.f_vs[i], x=self.f_vs[i]) / (self.ab_flux[i]/self.ab_source_flux)
            
            # If distance dependent
            if units in {'jy', 'fv', 'flux', 'fl'}:
                mag[i] = mag[i] / (4.0 * np.pi * d_L**2.0) 
            elif units == 'luminosity':
                mag[i] = mag[i]
            elif units == 'lsun':
                mag[i] = mag[i] / L_sun
            elif units == 'magnitude': 
                mag[i] = mag[i] / (4.0 * np.pi * d_L**2.0) / self.ab_source_flux   # normalized Flux        # Calculate the spectral flux density Fv in erg/s/cm^2/Hz
                mag[i] = -2.5 * np.log10(mag[i]) - 48.60 + app + to_vega + to_solar               
        
            else:
                raise NameError('Units of %s are unrecognized!' % units)

        if outspec_res is None: 
            self.spectrum['vs'] = vsn.reshape(1, vsn.size)
            self.spectrum['sed'] = sedn.T  # Changed to [particle, wavelength] format, need to add noise here!! 
        else:
            if isinstance(outspec_res, type(0.1)): 
                new_vsn = np.arange(vsn.min(), vsn.max(), (vsn.max()-vsn.min())/vsn.size/outspec_res)
            elif isinstance(outspec_res, type(np.array([0]))):
                new_vsn = outspec_res
            else:
                raise ValueError("The input outspec_res can only accept float value between 0 and 1 and a numpy arrary in Hertz to cover the interested spectrum energy range.", outspec_res)

            self.spectrum['vs'] = new_vsn.reshape(1, new_vsn.size)
            self.spectrum['sed'] = interp(new_vsn)  # Changed to [particle, wavelength] format, need to add noise here!! 

        return mag
