import os
from pymgal import utils
import re
from scipy.interpolate import interp1d
# import weight
# import collections
# import astro_filter_light
import numpy as np
import pyfits
from multiprocessing import Process, cpu_count, Queue, freeze_support


def worker(input, output):
    for func, args in iter(input.get, 'STOP'):
        result = calculate(func, args)
        output.put(result)


def calculate(func, args):
    return func(*args)


def fi(ages, seds, sage, smas):
    f = interp1d(ages, seds, bounds_error=False, fill_value="extrapolate")
    return f(sage) * smas


class SSP_models(object):
    r""" Load simple stellar population models.
    model = SSP_model(model_file, IMF="chab", metal=[], is_ised=False, is_fits=False,
    is_ascii=False, has_masses=False, units='a', age_units='gyrs', nsample=None)

    Parameters
    ----------
    model_file : File name for the SSP models. Only "name" + _ + "model",
                 such as c09_exp, or bc03_ssp. The whole name of the model
                 must be name_model_z_XXX_IMF.model (or .ised, .txt, .fits)

    IMF        : Name of the Initial mass function. Default: "chab"
                 Specify one of these IMF names: "chab", "krou" or "salp"

    metal      : The metallicity of the model. Specify at here if you only
                 want one or several metallicity included. Default: [],
                 all metallicity models are included to do interpolation.

    is_ised    : Is the model file an ised file? Default: False
    is_fits    : Is the model file a fits file? Default: False
    is_ascii   : Is the model file a ascii file? Default: False
    has_masses : Does it have mass? Default: False
    units      : The units of wavelength. Default: 'a'
    age_units  : The units of age. Default: 'gyrs'
                Note, this is the age units in the model file. It is changed into yrs in side the code.

    nsample    : The frequency sample points. Default: None, uses the model
                 frequences. Otherwise, the interpolated frequency will be used.
                 For the popurse of reducing required memory.
                 Int number or list, numpy array. If int number, frequency will be
                 assigned equally between min and max of model frequency.
                 If list, or numpy array, it will be directly interpolated with
                 the model frequency.
                 Note, must be consistent with the units.

    Returns
    -------
    Class for converting galaxy seds as a function of age to magnitudes
    as a function of redshift.

    See also
    --------
    get_seds function in this class

    Notes
    -----
    Reads in bc03 ised files as well as ascii files, effectively
    replacing cm_evolution from the bc03 src files

    Specify a model file. It should be a bc03 ised file, an ascii file,
    or a model file created with this class (which are stored as fits).
    See ezsps._load_ascii for information on formatting/units for ascii files
    This will automatically try to determine which type of file you have passed it.
    If this fails, you can specify it manually with the is_* flags.

    units, age_units, and has_masses apply only for ascii files.
    See ezsps._load_ascii()

    Example
    -------
    mm=pymgal.SSP_models("bc03_ssp", metal=[0.008])
    """
    tol = 1e-8  # tolerance for determining whether a given zf/z matches a stored zf/z

    def __init__(self, model_file='', IMF="chab", metal=[], is_ised=False,
                 is_fits=False, is_ascii=False, has_masses=False, units='a',
                 age_units='gyrs', nsample=None):

        self.nmets = len(metal)
        if self.nmets == 0:
            self.metals = np.array([])
        else:
            self.metals = np.array(metal)
        self.met_name = []

        self.zfs = None        # formation redshifts for models
        self.nzfs = []         # number of formation redshifts to model
        self.zs = None         # redshifts at which to project models
        self.nzs = []          # number of redshifts models should calculated at
        self.nfilters = []
        self.filters = None

        self.nages = []        # number of different aged SEDs
        self.ages = {}       # array of ages (years)
        # number of frequences in each SED array of frequences for SEDs number
        # of wavelengths in each SED (same as self.nvs)
        self.nvs = []
        self.vs = {}
        self.nls = []

        # array of wavelengths for SEDs (angstroms) age/SED grid. Units need
        # to be ergs/cm^2/s. Size is nvs x nages
        self.ls = {}
        self.seds = {}

        self.masses = {}
        self.sft = {}
        self.defined_freq = nsample

        self.model_dir = False
        if 'SSP_models' in os.environ:
            self.model_dir = os.environ['SSP_models']
        elif 'SSP_MODELS' in os.environ:
            self.model_dir = os.environ['SSP_MODELS']
        else:
            self.model_dir = os.path.dirname(
                os.path.realpath(__file__)) + '/models/'
            # save path to data folder: module directory/data

        if self.model_dir and (self.model_dir[-1] != os.sep):
            self.model_dir += os.sep

        # load model
        self._load(model_file, IMF, is_ised=is_ised, is_fits=is_fits,
                   is_ascii=is_ascii, has_masses=has_masses, units=units,
                   age_units=age_units)

    def _load(self, model_file, IMF, is_ised=False, is_fits=False,
              is_ascii=False, has_masses=False, units='a', age_units='gyrs'):

        # find the model file
        self._find_model_file(model_file, IMF)

        model_file = self.filenames[0]

        # Load input file depending on what type of file it is.
        # test for a bruzual-charlot binary ised file
        if model_file[len(model_file) - 5:] == '.ised' or is_ised:
            self._load_ised(self.filenames)
        else:
            # And then the rest.
            # fp = open(model_file, 'r')
            # start = fp.read(80)
            # fp.close()
            start = np.fromfile(model_file, dtype="c", count=80).tostring()
            if (re.search('SIMPLE\s+=\s+T', start.decode('utf-8'), re.IGNORECASE) or is_fits) and not(is_ascii):
                self._load_model(self.filenames)
            else:
                self._load_ascii(self.filenames, has_masses=has_masses,
                                 units=units, age_units=age_units)

        if self.defined_freq is not None:
            if isinstance(self.defined_freq, type(1)):
                for i, mn in enumerate(self.met_name):
                    newfreq = np.arange(self.vs[mn].min(), self.vs[mn].max(
                    ), (self.vs[mn].max() - self.vs[mn].min() * 0.999) / self.defined_freq)
                    f = interp1d(self.vs[mn], self.seds[mn].T,
                                 bounds_error=False, fill_value="extrapolate")
                    self.seds[mn] = f(newfreq).T
                    self.vs[mn] = newfreq
                    self.nvs[i] = self.nls[i] = self.defined_freq
                    self.ls[mn] = utils.to_lambda(newfreq)
            elif isinstance(self.defined_freq, type([])) or isinstance(self.defined_freq, type(np.ones(1))):
                newfreq = np.asarray(self.defined_freq)
                for i, mn in enumerate(self.met_name):
                    f = interp1d(self.vs[mn], self.seds[mn].T,
                                 bounds_error=False, fill_value="extrapolate")
                    self.seds[mn] = f(newfreq).T
                    self.vs[mn] = newfreq
                    self.nvs[i] = self.nls[i] = len(newfreq)
                    self.ls[mn] = utils.to_lambda(newfreq)
            else:
                raise ValueError(
                    "Do not know what is this input frequency %s " % self.defined_freq)

        if has_masses:  # need to normalize the seds.
            for metmodel in self.met_name:
                # because the fist mass is zero...
                self.seds[metmodel][:, 1:] /= self.masses[metmodel][1:]

        # always include a t=0 SED to avoid out of age interpolation errors later
        # a t=0 SED is also assumed during CSP generation
        # if self.nages and self.ages.min() > 0:
        #     self.ages = np.append(0, self.ages)
        #     self.seds = np.hstack((np.zeros((self.nvs, 1)), self.seds))
        #     if self.has_masses:
        #         self.masses = np.append(0, self.masses)
        #     self.nages += 1

    def _find_model_file(self, file, IMF):
        # first check file path
        files = []
        metals = []

        for filen in os.listdir(self.model_dir):
            if (filen[:len(file)] == file) and (filen[-10:-6] == IMF):
                filemet = np.float64(filen.split("_")[-2])
                if self.nmets == 0:
                    files.append('%s%s' % (self.model_dir, filen))
                    metals.append(filemet)
                    self.met_name.append(filen.split("_")[-2])
                else:
                    if np.min(np.abs(self.metals - filemet)) <= self.tol * 10**4:
                        files.append('%s%s' % (self.model_dir, filen))
                        metals.append(filemet)
                        self.met_name.append(filen.split("_")[-2])

        # now loop through the different files
        if len(files) != 0:
            for check in files:
                if os.path.isfile(check):
                    print("Searching SSP model files, found: %s" % check)
                else:
                    print("There is something wrong with this model file: %s", check)
                    print("!! This may cause problems in reading !!")
            self.filenames = files
            # force to use model metals
            self.metals = np.array(metals, dtype=np.float64)
            self.nmets = self.metals.size
        else:
            raise ValueError(
                'There is no such file/files: %s was not found!' % file)

    def _load_ised(self, file):
        """
        ezsps._load_ised(file)

        Load a bruzual and charlot binary ised file.
        Saves the model information in the model object
        """

        for i in range(self.nmets):
            (seds, ages, vs) = utils.read_ised(file[i])
            # store ages
            self.ages[self.met_name[i]] = ages
            self.nages.append(ages.size)
            # store wavelength
            self.vs[self.met_name[i]] = vs
            self.nvs.append(vs.size)
            self.ls[self.met_name[i]] = utils.to_lambda(self.vs)
            self.nls.append(vs.size)
            # store seds
            self.seds[self.met_name[i]] = seds

    def _load_model(self, model_file):
        """
        ezsps._load_model(model_file)

        loads a model from a fits file created with ezsps.save_model()
        Saves the model information in the model object
        """

        for i in range(self.nmets):
            fits = pyfits.open(model_file[i])
            # was sed information included in this model file?
            if fits[0].header['has_seds']:
                # store seds
                self.seds[self.met_name[i]] = fits[0].data
                # store wavelength
                self.vs[self.met_name[i]] = fits[1].data.field('vs')
                self.nvs.append(fits[1].data.field('vs').size)
                self.ls[self.met_name[i]] = utils.to_lambda(fits[1].data.field('vs'))
                self.nls.append(fits[1].data.field('vs').size)
                # store ages
                self.ages[self.met_name[i]] = fits[2].data.field('ages')
                self.nages.append(fits[2].data.field('ages').size)
                # how about masses?
                if 'has_mass' in fits[2].header and fits[2].header['has_mass']:
                    self.has_masses = True
                    self.masses[self.met_name[i]] = fits[2].data.field('masses')
                # and sfh?
                if 'has_sfh' in fits[2].header and fits[2].header['has_sfh']:
                    self.sfh[self.met_name[i]] = fits[2].data.field('sfh')
                    self.has_sfh = True
            else:
                raise ValueError(
                    'Fits file: %s does not contains seds!' % file)

    def _load_ascii(self, file, has_masses=False, units='a', age_units='gyrs'):
        """
        ezsps._load_ascii(file, has_masses=False, units='a', age_units='gyrs')

        Load a model file in ascii format.
        The file should be a data array of size (nwavelengths+1,nages+1).
        The first row specifies the age of each column, and the first column
        specifies the wavelength of each row.
        This means the data value in the first row of the first column is ignored.
        However, it still must have SOME value ('0' is fine) as a placeholder
        You can include masses in the file by specifying the mass (in solar masses)
        for each age in the second row.
        If you do this then set has_masses=True
        It loads a bit slow, so you should save it as a fits -
        see ezsps.save_model() -
        if you are going to be using it more than once.

        Specify units for the age with 'age_units'.  Default is gyrs.
        See ezsps.utils.to_years() for avaialable unit specifications
        Specify units for wavelength & flux with 'units'.
        Default is 'a' for Angstroms, with flux units of ergs/s/angstrom.
        Set units='hz' for frequency with flux units of 'ergs/s/hertz'
        You can also set the units as anything else found in ezsps.utils.to_meters()
        as long as the flux has units of ergs/s/(wavelength units)
        """

        for i in range(self.nmets):
            model = utils.rascii(file[i])

            self.vs[self.met_name[i]] = model[1:, 0]
            self.nvs.append(model[1:, 0].size)
            self.nls.append(model[1:, 0].size)
            self.ages[self.met_name[i]] = model[0, 1:]
            self.nages.append(model[0, 1:].size)

            if has_masses:
                self.masses[self.met_name[i]] = model[1, 1:]
                self.seds[self.met_name[i]] = model[1:, 2:]
            else:
                self.seds[self.met_name[i]] = model[1:, 1:]

            # convert to intermediate units (hz, ergs/s/hz)
            units = units.lower()
            age_units = age_units.lower()
            if units != 'hz':
                self.seds[self.met_name[i]] *= self.vs[self.met_name[i]].reshape(
                    (self.nvs[-1], 1))**2.0 / utils.convert_length(utils.c, outgoing=units)
                self.vs[self.met_name[i]] = utils.to_hertz(
                    self.vs[self.met_name[i]], units=units)

            self.ls[self.met_name[i]] = utils.to_lambda(
                self.vs[self.met_name[i]])

            # convert from ergs/s/Hz to ergs/s/Hz/cm^2.0 @ 10pc
            self.seds[self.met_name[i]] /= 4.0 * np.pi * \
                utils.convert_length(10, incoming='pc', outgoing='cm')**2.0

            # convert ages to the proper units
            self.ages[self.met_name[i]] = utils.to_years(
                self.ages[self.met_name[i]], units=age_units)

            # now sort it to make sure that age is increasing
            # sind = self.ages[self.met_name[i]].argsort()
            # self.ages = self.ages[sind]
            # self.seds = self.seds[:, sind]

            # # the same for frequency
            # sind = self.vs.argsort()
            # self.vs = self.vs[sind]
            # self.seds = self.seds[sind, :]

    def get_seds(self, simdata, Ncpu=None, dust_func=None, units='Fv'):
        r"""
        Seds = SSP_model(simdata, Ncpu=None, dust_func=None, units='Fv')

        Parameters
        ----------
        simudata   : Simulation data read from load_data
        Ncpu       : The number of CPUs for parallel interpolation
        dust_func  : dust function.
        units      : The units for retruned SEDS. Default: 'Fv'

        Get SEDS for the simulated star particles

        Available output units are (case insensitive):

        ========== ====================
        Name       Units
        ========== ====================
        Jy         Jansky
        Fv         ergs/s/cm^2/Hz
        Fl         ergs/s/cm^2/Angstrom
        Flux       ergs/s/cm^2
        Luminosity ergs/s
        ========== ====================
        """

        seds = np.zeros((self.nvs[0], simdata.S_age.size), dtype=np.float64)

        # We do not do 2D interpolation since there is only several metallicity
        # in the models.
        if self.nmets > 1:
            mids = np.interp(simdata.S_metal, self.metals, np.arange(self.nmets))
            mids = np.int32(np.round(mids))

        for i, metmodel in enumerate(self.met_name):
            # f = interp1d(self.ages[metmodel], self.seds[metmodel],
            #              bounds_error=False, fill_value="extrapolate")

            if self.nmets > 1:
                ids = np.where(mids == i)[0]
            else:
                ids = np.ones(simdata.S_metal.size, dtype=np.bool)

            # parallel
            freeze_support()
            if Ncpu is None:
                NUMBER_OF_PROCESSES = cpu_count()
            else:
                NUMBER_OF_PROCESSES = Ncpu

            Ns = 20000
            # Number of total Tasks, control the size of passing arrays
            N = np.int32(ids.size / Ns)
            if N < NUMBER_OF_PROCESSES:
                N = NUMBER_OF_PROCESSES
                Ns = np.int32(ids.size / N)

            # Create queues
            task_queue = Queue()
            done_queue = Queue()

            Tasks = [(fi, (self.ages[metmodel], self.seds[metmodel],
                           simdata.S_age[ids][i * Ns:(i + 1) * Ns],
                           simdata.S_mass[ids][i * Ns:(i + 1) * Ns])) for i in range(N)]
            if ids.size - N * Ns > 0:
                Tasks.append((fi, (self.ages[metmodel], self.seds[metmodel],
                                   simdata.S_age[ids][N * Ns:ids.size],
                                   simdata.S_mass[ids][N * Ns:ids.size])))

            # Submit tasks
            for task in Tasks:
                task_queue.put(task)

            # Start worker processes
            for i in range(NUMBER_OF_PROCESSES):
                Process(target=worker, args=(task_queue, done_queue)).start()

            # Get results
            for i in range(len(Tasks)):
                if i < N:
                    seds[:, ids[i * Ns:(i + 1) * Ns]] = done_queue.get()
                elif i == N:
                    seds[:, ids[N * Ns:]] = done_queue.get()
                else:
                    raise ValueError("Wrong in number of Tasks %d, %d" % (i, len(Tasks)))

            # Tell child processes to stop
            for i in range(NUMBER_OF_PROCESSES):
                task_queue.put('STOP')

            if dust_func is not None:
                seds[:, ids] *= dust_func(simdata.S_age[ids], self.ls[metmodel])
            #     seds[:, ids] = f(simdata.S_age[ids]) * simdata.S_mass[ids]
            # else:
                # seds[:, ids] = f(simdata.S_age[ids]) * simdata.S_mass[ids] * \
                #     dust_func(simdata.S_age[ids], self.ls[metmodel])

        # if simdata.grid_mass is not None:
        #     seds = binned_statistic_2d(simdata.S_pos[:, 0], simdata.S_pos[:, 1],
        #                                values=seds,
        #                                bins=[simdata.nx, simdata.nx],
        #                                statistic='sum')[0]

        units = units.lower()
        if units == 'jy':
            return seds / 1e-23
        if units == 'fv':
            return seds
        if units == 'fl':
            seds *= (self.vs[self.met_name[0]].reshape(self.nvs[0], 1)
                     )**2.0 / utils.convert_length(utils.c, outgoing='a')
            return seds

        seds *= self.vs[self.met_name[0]].reshape(self.nvs[0], 1)
        if units == 'flux':
            return seds
        if units == 'luminosity':
            return seds * 4.0 * np.pi * utils.convert_length(10, incoming='pc', outgoing='cm')**2.0
        raise NameError('Units of %s are unrecognized!' % units)
