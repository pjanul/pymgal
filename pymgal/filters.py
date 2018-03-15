import os
import numpy as np
import pyfits
from pymgal import utils
from scipy.integrate import simps
from scipy.interpolate import interp1d
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
        else:
            self.f_vs[fname] = utils.to_hertz(filters[:, 0], units=units)
            self.f_ls[fname] = utils.convert_length(
                filters[:, 0], incoming=units, outgoing='a')
        self.f_tran[fname] = filters[:, 1]
        self.npts[fname] = filters[:, 0].size

        # normalization for calculating ab mags for this filter
        self.ab_flux[fname] = self.ab_source_flux * \
            simps(self.f_tran[fname] / self.f_vs[fname], self.f_vs[fname])

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
        sed_flux = simps(interp(self.f_vs[fn]).T * self.f_tran[fn] / self.f_vs[fn], self.f_vs[fn])
        return -2.5 * np.log10(sed_flux / self.ab_flux[fn])

    ##############
    #  calc mag  #
    ##############
    def calc_mag(self, sspmod, simd, fn=None, Ncpu=None, dust_func=None, redshift=None,
                 apparent=False, vega=False, flux=False):
        r"""
        mag = pymgal.calc_mag(ssp, simd, z, fn=None)

        sspmod      : loaded ssp models.
        simd        : loaded simulation data from load_data
        fn          : the name of filter/s. string, list of strings, or None
                        By default (None), all loaded filters will be included in the calculation
        Ncpu        : The number of CPUs for parallel interpolation
        dust_func   : The function for dust attenuetion.
        redshift    : Specified redshift of the object. Default: None, the simulation redshift will
                        be used.
        apparent    : If you need apparent magnitude, set True. Default False
        vega        : AB magnitude in VEGA.
        flux        : return flux.
        returns     :   0. Absolute AB magnitude Default
                        1. apparent magnitude if apparent=True
                        2. Absolute AB magnitude Outputs vega mags if vega=True
                        3. flux, if flux=True

        Example:
            mag = pymgal.calc_mag(ssp, simd, z, fn=None)

        Notes:
        Calculate the absolute AB magnitude (or VEGA) of the given sed at the given redshift.
        Set ``z=0`` for rest-frame magnitudes.
        ``vs`` should give the frequency (in Hz) of every point in the SED,
        and the sed should have units of ergs/s/cm^2/Hz.
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

        if not redshift:
            redshift = simd.z

        if apparent:
            app = 5. * np.log10(simd.cosmology.luminosity_distance(redshift).value
                                * 1.0e5) if redshift > 0 else -np.inf
        else:
            app = 0.0

        vs = sspmod.vs[sspmod.met_name[0]]
        mag = {}
        shifted = vs / (1 + redshift)

        interp = interp1d(vs, sspmod.get_seds(simd, Ncpu=Ncpu, dust_func=dust_func), axis=0,
                          bounds_error=False, fill_value="extrapolate")

        for i in fn:
            if vega:
                to_vega = self.vega_mag[i]
            else:
                to_vega = 0.0

            c = ((shifted > self.f_vs[i].min()) & (shifted < self.f_vs[i].max())).sum()
            if c < 3:
                print("Warning, wavelength range from SSP models is too near filter,",
                      " magnitude is assigned nan")
                mag[i] = np.nan
            # and that the SED actually covers the whole filter
            if shifted.min() > self.f_vs[i].min() or shifted.max() < self.f_vs[i].max():
                print("Warning, wavelength range from SSP models is outside of filter,",
                      " magnitude is assigned nan")
                mag[i] = np.nan

            if flux:
                mag[i] = (1 + redshift) * simps(interp(self.f_vs[i] * (1 + redshift)).T *
                                                self.f_tran[i] / self.f_vs[i], self.f_vs[i]) / self.ab_flux[i]

            else:
                sed_flux = (1 + redshift) * simps(interp(self.f_vs[i] * (1 + redshift)).T *
                                                  self.f_tran[i] / self.f_vs[i], self.f_vs[i])
                mag[i] = -2.5 * np.log10(sed_flux / self.ab_flux[i]) + app + to_vega
        return mag
