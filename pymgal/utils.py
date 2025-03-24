import os
import re
import numpy as np
from numba import njit, prange
from scipy.spatial import KDTree

# some useful constants
c = 299792458            # speed of light (m/sec)
m_per_au = 1.49598e11        # meters per astronomical unit
au_per_pc = 3600 * 180 / np.pi    # AUs per parsec


def to_years(to_convert, units='gyrs', reverse=False):
    """
    res = ezsps.utils.to_years(to_convert, units='gyrs', reverse=False)

    Converts the given time to years from the given units.
    If reverse=True then it converts from years into the given units

    :param to_convert: The time to convert
    :param units: The units to convert the time to
    :param reverse: Converts from years if True
    :type to_convert: int, float
    :type units: string
    :type reverse: bool
    :returns: The converted time
    :rtype: int, float

    :Example:
        >>> import ezsps
        >>> ezsps.utils.to_years(1e-9, units='gyrs')
        1.0

    **units**

    Available units are (case insensitive):

    ==================  ====================
           Name             Units
    ==================  ====================
    gigayears,gyrs,gyr  Gigayears
    megayears,myrs,myr  Megayears
    years,yrs,yr        Years
    days,day            Days
    seconds,secs,sec,s  Seconds
    log                 log10(years)
    ==================  ====================

    .. seealso:: :func:`ezsps.utils.convert_time`
    """

    units = units.lower()
    factor = 0
    to_convert = np.asarray(to_convert)
    if units == 'gigayears' or units == 'gyrs' or units == 'gyr':
        factor = 1e9
    if units == 'megayears' or units == 'myrs' or units == 'myr':
        factor = 1e6
    if units == 'years' or units == 'yrs' or units == 'yr':
        factor = 1
    if units == 'days' or units == 'day':
        factor = 1.0 / 365.0
    if units == 'seconds' or units == 'secs' or units == 'sec' or units == 's':
        factor = 1.0 / (365.0 * 86400)
    if factor != 0:
        if reverse:
            factor = 1.0 / factor
        return to_convert * factor

    if units == 'log':
        if reverse:
            return np.log10(to_convert)
        return 10.0**to_convert

    raise NameError('Units of %s are not recognized!' % units)


def convert_time(to_convert, incoming='secs', outgoing='gyrs'):
    """
    res = ezsps.utils.convert_time(to_convert, incoming='secs', outgoing='gyrs')

    Converts the given time from the incoming units to the outgoing units.

    :param to_convert: The length to convert
    :param incoming: The units to convert the time from
    :param outgoing: The units to convert the time to
    :type to_convert: int, float
    :type incoming: string
    :type outgoing: string
    :returns: The converted time
    :rtype: int, float

    :Example:
        >>> import ezsps
        >>> ezsps.utils.convert_time(1, incoming='years', outgoing='s')
        31536000.0

    .. seealso:: see :func:`ezsps.utils.to_years` for available units."""

    return to_years(to_years(to_convert, units=incoming), units=outgoing, reverse=True)


def to_meters(to_convert, units='a'):
    """
    res = ezsps.utils.to_meters(to_convert, units='a')

    Converts a length from the given units to meters

    :param to_convert: The length to convert
    :param units: The units to convert the length to
    :type to_convert: int, float
    :type units: string
    :returns: The converted length
    :rtype: int, float

    :Example:
        >>> import ezsps
        >>> ezsps.utils.to_meters(1e10, units='a')
        1.0

    **units**

    Available units are (case insensitive):

    ================= ====================
           Name             Units
    ================= ====================
    a,angstroms       Angstroms
    nm,nanometers     Nanometers
    um,microns        Microns
    mm,millimeters    Millimeters
    cm,centimeters    Centimeters
    m,meters          Meters
    km,kilometers     Kilometers
    au                Astronomical Units
    pc,parsecs        Parsecs
    kpc, kiloparsecs  Kiloparsecs
    mpc, megaparsecs  Megaparsecs
    ================= ====================

    .. seealso:: :func:`ezsps.utils.convert_length`

    """

    units = units.lower()
    to_convert = np.asarray(to_convert)
    if units == 'angstroms' or units == 'a':
        return to_convert / 1e10
    if units == 'nanometers' or units == 'nm':
        return to_convert / 1e9
    if units == 'microns' or units == 'um':
        return to_convert / 1e6
    if units == 'millimeters' or units == 'mm':
        return to_convert / 1e3
    if units == 'centimeters' or units == 'cm':
        return to_convert / 1e2
    if units == 'meters' or units == 'm':
        return to_convert
    if units == 'kilometers' or units == 'km':
        return to_convert * 1000.0
    if units == 'au':
        return to_convert * m_per_au
    if units == 'parsecs' or units == 'pc':
        return to_convert * m_per_au * au_per_pc
    if units == 'kilparsecs' or units == 'kpc':
        return to_convert * m_per_au * au_per_pc * 1000.0
    if units == 'megaparsecs' or units == 'mpc':
        return to_convert * m_per_au * au_per_pc * 1e6

    raise NameError('Units of %s are not recognized!' % units)


def to_hertz(to_convert, units='a'):
    """
    res = ezsps.utils.to_hertz(to_convert, units='Angstroms')

    Converts the given wavelength (in the given units) to hertz.

    :param to_convert: The wavelength to convert
    :param units: The units the wavelength is in
    :type to_convert: int, float
    :type units: string
    :returns: The converted frequency
    :rtype: float

    :Example:
        >>> import ezsps
        >>> ezsps.utils.to_hertz(1000, units='a')
        2997924580000000.0

    .. seealso::
        see :func:`ezsps.utils.to_meters` for list of available units

        Also see :func:`ezsps.utils.to_lambda`
    """

    return (c / to_meters(1.0, units=units)) / np.asarray(to_convert)


def to_lambda(to_convert, units='a'):
    """
    res = ezsps.utils.to_lambda(to_convert, units='a')

    Converts the given frequency to a wavelength in the given units.

    :param to_convert: The frequency to convert
    :param units: The desired units of the output wavelength
    :type to_convert: int, float
    :type units: string
    :returns: The converted wavelength
    :rtype: float

    :Example:
        >>> import ezsps
        >>> ezsps.utils.to_lambda(2997924580000000.0, units='a')
        1000.0

    .. seealso::
        see :func:`ezsps.utils.to_meters` for list of available units

        Also see :func:`ezsps.utils.to_hertz`
    """

    return (c / to_meters(1.0, units=units)) / np.asarray(to_convert)


def convert_length(to_convert, incoming='m', outgoing='a'):
    """
    res = ezsps.utils.convert_length(to_convert, incoming='m', outgoing='a')

    converts a length from the incoming units to the outgoing units.

    :param to_convert: The length to convert
    :param incoming: The units to convert the length from
    :param outgoing: The units to convert the length to
    :type to_convert: int, float
    :type incoming: string
    :type outgoing: string
    :returns: The converted length
    :rtype: int, float

    :Example:
        >>> import ezsps
        >>> ezsps.utils.convert_length(1, incoming='pc', outgoing='au')
        206264.80624709636

    .. seealso:: see :func:`ezsps.utils.to_meters` for available units.
    """

    return to_meters(to_convert, units=incoming) / to_meters(1.0, units=outgoing)


def rascii(filename, silent=False):
    """
    res = ezsps.utils.rascii(filename, silent=False)

    Reads in numeric data stored in an ascii file into a numpy array.

    :param filename: The name of the ascii file
    :param silent: Whether or not to output basic file information
    :type filename: string
    :type silent: bool
    :returns: A numpy array
    :rtype: np.array()

    .. warning::
        Skips any lines that have any non-numeric data, and any data
        lines with a different number of columns than the first data line.

    .. seealso:: :func:`ezsps.utils.wascii`
    """
    file = open(filename, 'r')
    found = False
    nlines = 0
    ngood = 0

    for line in file:

        nlines += 1
        if re.search('^\s*$', line) or re.search('[^\s\d.eEdD\-+]', line):
            continue

        parts = line.split()
        nparts = len(parts)

        if not found:
            found = True
            allowed = nparts
            res = parts
            continue

        if nparts != allowed:
            continue

        ngood += 1
        res.extend(parts)

    ngood += 1
    if ngood == 0:
        return np.array([])

    arr = np.array(res)
    arr.shape = (ngood, -1)

    return np.array(res).reshape(ngood, -1).astype('float')


def wascii(array, filename, formats, blank=False, header=None, names=None):
    """
    ezsps.utils.wascii(array, filename, formats, blank=False, header=False, names=None)

    Writes out a np array to a well formated file.

    :param array: The numpy array to write out
    :param filename: The name of the output file
    :param formats: A list of python string formats (one for each column)
    :param blank: Whether or not to output a blank line at the end of the file
    :param header: A string or list of strings to write out as the header
    :param names: A list of column names with which to build a header
    :type array: a 2D numpy array
    :type filename: string
    :type formats: string,list
    :type blank: bool
    :type header: string,list
    :type blank: string,list
     """

    table = np.asarray(array)
    if table.ndim != 2:
        raise NameError('I was expecting a 2D data table')
    nrows, ncols = table.shape

    if type(formats) is str:
        formats = [formats] * ncols

    if ncols != len(formats):
        raise NameError(
            'Number of supplied formats does not match number of table columns!')

    # if column names were provided, create a header that list column
    # names/numbers
    if names is not None:
        if len(names) != ncols:
            raise NameError(
                'Number of supplied column names does not match number of table columns!')
        if header is None:
            header = []
        header.append('# Column Descriptions:')
        name_format = '# %0' + \
            ('%1d' % (np.ceil(np.log10(ncols)))) + 'd: %s'
        for i in range(ncols):
            header.append(name_format % (i + 1, names[i]))

    if (header is not None) & (not isinstance(header, type(''))):
        header = "\n".join(header)

    if ncols == 1:
        file = "\n".join(formats[0] % val for val in table.ravel())
    else:
        strings = [''] * nrows
        for i in range(nrows):
            strings[i] = ' '.join(
                [format % val for format, val in zip(formats, table[i, :])])
        file = "\n".join(strings)

    fh = open(filename, 'wb')
    if header is not None:
        fh.write(header + "\n")
    fh.write(file)
    if blank:
        fh.write("\n")
    fh.close()


def _read_binary(fhandle, type='i', number=1, swap=False):
    """
    res = ezsps.utils._read_binary(fhandle, type='i', number=1, swap=False)

    reads 'number' binary characters of type 'type' from file handle 'fhandle'
    returns the value (for one character read) or a numpy array
    set swap=True to byte swap the array after reading
    """

    return np.fromfile(fhandle, dtype=type, count=number)

    # arr = array.array(type)
    # arr.fromfile(fhandle, number)
    # if swap:
    #     arr.byteswap()
    # if len(arr) == 1:
    #     return arr[0]
    # else:
    #     return np.asarray(arr)


def read_ised(file):
    """ (seds, ages, vs) = ezsps.utils.read_ised(file)

    Read a bruzual and charlot binary ised file.

    :param file: The name of the ised file
    :type file: string
    :returns: A tuple containing model data
    :rtype: tuple

    .. note::
        All returned variables are numpy arrays.
        ages and vs are one dimensional arrays,
        and seds has a shape of (vs.size,ages.size)

    **units**
    Returns units of:

    =============== ===============
    Return Variable   Units
    =============== ===============
    seds            Ergs/s/cm**2/Hz
    ages            Years
    vs              Hz
    =============== ===============
    """

    if not(os.path.isfile(file)):
        raise ValueError('The specified model file was not found!')

    # open the ised file
    fh = open(file, 'rb')

    # start reading
    junk = _read_binary(fh)
    nages = _read_binary(fh)

    # first consistency check
    if nages < 1 or nages > 2000:
        raise ValueError(
            'Problem reading ised file - unexpected data found for the number of ages!')

    # read ages
    ages = np.asarray(_read_binary(fh, type='f', number=nages))

    # read in a bunch of stuff that I'm not interested in but which I read
    # like this to make sure I get to the right spot in the file
    junk = _read_binary(fh, number=2)
    iseg = _read_binary(fh, number=1)
    if iseg > 0:
        junk = _read_binary(fh, type='f', number=6 * iseg)
    junk = _read_binary(fh, type='f', number=3)
    junk = _read_binary(fh)
    junk = _read_binary(fh, type='f')
    junk = _read_binary(fh, type='c', number=80)
    junk = _read_binary(fh, type='f', number=4)
    junk = _read_binary(fh, type='c', number=160)
    junk = _read_binary(fh)
    junk = _read_binary(fh, number=3)

    # read in the wavelength data
    nvs = _read_binary(fh)

    # consistency check
    if nvs < 10 or nvs > 12000:
        raise ValueError(
            'Problem reading ised file - unexpected data found for the number of wavelengths!')

    # read wavelengths and convert to frequency (comes in as Angstroms)
    # also reverse the array so it will be sorted after converting to frequency
    ls = _read_binary(fh, type='f', number=nvs)[::-1]

    # create an array for storing SED info
    seds = np.zeros((nvs, nages))

    # now loop through and read in all the ages
    for i in range(nages):
        junk = _read_binary(fh, number=2)
        nv = _read_binary(fh)
        if nv != nvs:
            raise ValueError(
                'Problem reading ised file - unexpected data found while reading seds!')

        seds[:, i] = _read_binary(fh, type='f', number=nvs)[::-1]
        nx = _read_binary(fh)
        junk = _read_binary(fh, type='f', number=nx)
        del(junk)
    # now convert the seds from Lo/A to ergs/s/Hz
    seds *= 3.826e33 * ls.reshape((nvs, 1))**2.0 / \
        convert_length(c, outgoing='a')
    # convert from ergs/s/Hz to ergs/s/Hz/cm^2.0 @ 10pc
    seds /= 4.0 * np.pi * convert_length(10, incoming='pc', outgoing='cm')**2.0
    vs = to_hertz(ls)

    fh.close()

    # sort in frequency space
    sinds = vs.argsort()

    return (seds[sinds, :], ages, vs[sinds, :])


def numba_interp1d(x, y):
    """
    Returns a function that performs linear interpolation with Numba acceleration,
    matching scipy's interp1d behavior with bounds_error=False and fill_value="extrapolate".
    
    Parameters:
      - x: 1D array of known x data (must be sorted).
      - y: 2D array of known y data (interpolation happens along axis=1).
    
    Returns:
      - A function that takes an array xnew and returns the interpolated y values.
    """
    @njit(parallel=True)
    def interp_function(xnew):
        npts = len(xnew)
        ynew = np.empty((y.shape[0], npts), dtype=y.dtype)
        for i in prange(npts):
            xi = xnew[i]
            # Extrapolation: use first interval if xi is below the known range
            if xi < x[0]:
                j = 0
                t = (xi - x[0]) / (x[1] - x[0])
            # Extrapolation: use last interval if xi is above the known range
            elif xi > x[-1]:
                j = len(x) - 2
                t = (xi - x[-2]) / (x[-1] - x[-2])
            else:
                # xi is within the known range.
                # Special case: if xi exactly equals the last x value, then return y[:, -1]
                if xi == x[-1]:
                    j = len(x) - 2
                    t = 1.0
                else:
                    # Find the interval such that x[j] <= xi < x[j+1]
                    for j_temp in range(len(x) - 1):
                        if xi >= x[j_temp] and xi < x[j_temp+1]:
                            j = j_temp
                            t = (xi - x[j]) / (x[j+1] - x[j])
                            break

            # Compute linear interpolation (or extrapolation)
            for k in range(y.shape[0]):
                ynew[k, i] = (1 - t) * y[k, j] + t * y[k, j+1]
        return ynew
    return interp_function


@njit
def scale_seds(seds, S_mass, ids):
    seds[:, ids] *= S_mass[ids]
    return seds

@njit
def apply_dust(seds, ids, dust_arr):
    seds[:, ids] *= dust_arr
    return seds

# Create kD tree and extract the distance of the kth nearest neighbour (kth distance = 1 sigma of gaussian kernel)
def knn_distance(positions, k):
    kdtree = KDTree(positions)
    _, indices = kdtree.query(positions, k + 1)  # k + 1 to exclude the point itself
    kth_neighbors = np.take(positions, indices[:, k], axis=0)
    kth_distances = np.linalg.norm(kth_neighbors - positions, axis=1)
    return kth_distances


