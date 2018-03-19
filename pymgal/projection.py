import numpy as np
import astropy.units as u
import re
from astropy.coordinates import SkyCoord
from astropy.time import Time
import version
# scipy must >= 0.17 to properly use this!
# from scipy.stats import binned_statistic_2d


def get_property(prop):
    result = re.search(r'{}\s*=\s*[\'"]([^\'"]*)[\'"]'.format(prop), open('__init__.py').read())
    return result.group(1)


class projection(object):
    r"""load analysing data from simulation snapshots (gadget format only),
    yt, or raw data. Currently only works with snapshot=True. Will work more
    on other data sets.

    Parameters
    ----------
    wdata   : The data to be saved. Type: dictionary of array.
                It is coming from the outputs of filters.calc_mag.
                The array must be the same length of simulation data or npx x npx.
    simd    : The loaded simulation data from load_data function.
    axis    : can be 'x', 'y', 'z', or a *list* of degrees [alpha, beta, gamma],
                which will rotate the data points by $\alpha$ around the x-axis,
                $\beta$ around the y-axis, and $\gamma$ around the z-axis.
                or direction of a vector pointing to the line of sight in np.array,
                or directly the rotation matrix [3x3] in np.array.
                Default: "z"
    npx     : The pixel number of the grid for projection.
                Type: int. Default: 512
                A npx x npx image will be produced later.
    AR      : Angular resolution. Type: arcsec. Default: None
    redshift: The redshift of the object at. Default: None.
                If None, redshift from simulation data will be used.
                Maybe move to 0.01 if it is 0.
    zthick  : The thickness in projection direction. Default: None.
                If None, use all data from cutting region. Otherwise set a value in simulation
                length unit (kpc/h normally), then a slice of data [center-zthick, center+zthick]
                will be used to make the y-map.
    SP      : Faked sky positions in [RA (longitude), DEC (latitude)] in degrees.
                Default: [194.95, 27.98], Coma' position.
                If [x,y,z] (len(SP) == 3) of the Earth position in the pos coordinate is given,
                The pos - [x,y,z] are taken as the J2000 3D coordinates and converted into RA, DEC.
    flux    : is the input wdata in luminosity/flux? Default: False, assume in ab mag.
                Set this to true if particles' luminosity are given in wdata.
    outmas  : do you want to out put stellar mass? Default: False.
                If True, the stellar mass in each pixel are saved.
    outage  : do you want to out put stellar age (mass weighted)? Default: False.
                If True, the stellar age in each pixel are saved.
    outmet  : do you want to out put stellar metallicity (mass weighted)? Default: False.
                If True, the stellar metallicity in each pixel are saved.
    Note    : The comments in str will be put into the fit file header. Defualt: None
    Notes
    -----


    Example
    -------
    Pdata = pymgal.projection(part_lum, simu_data, npx=1024, flux=Ture)
    Pdata.write_fits_image("filename.fits")
    """

    def __init__(self, data, simd, axis="z", npx=512, AR=None, redshift=None, zthick=None,
                 SP=[194.95, 27.98], flux=False, outmas=False, outage=False, outmet=False,
                 Note='None'):

        self.axis = axis
        self.npx = npx
        self.ar = AR
        if redshift is None:
            self.z = simd.z
        else:
            self.z = redshift
        self.flux = flux

        self.pxsize = 0.
        if len(SP) == 2:
            self.sp = SP
            self.cc = simd.center
        elif len(SP) == 3:
            self.sp = False
            self.cc = SP
        else:
            raise ValueError("SP length should be either 2 or 3!")
        self.rr = simd.radius
        self.flux = flux
        self.omas = outmas
        self.oage = outage
        self.omet = outmet
        self.zthick = zthick
        self.outd = {}
        self.note = Note

        # if not flux:
        #     if isinstance(data, type({})):
        #         self.outd = {}
        #         for i in data.keys():
        #             data[i] = 10**(data[i]/-2.5)
        #             # self.outd[i] = np.zeros((npx, npx), dtype=float)
        #     else:
        #         raise ValueError("Do not accept this data type %s " % type(data))

        self._prep_out(data, simd)

    def _prep_out(self, d, s):
        r""" rotate the data points and project them into a 2D grid.
        """
        # ratation data points first

        pos = np.copy(s.S_pos)
        pos -= self.cc

        if isinstance(self.axis, type('')):
            if self.axis.lower() == 'y':  # x-z plane
                pos = pos[:, [0, 2, 1]]
                if self.sp is False:
                    self.cc[0], self.cc[1], self.cc[2] = s.center[0] - self.cc[0],\
                        s.center[2] - self.cc[2], s.center[1] - self.cc[1]
            elif self.axis.lower() == 'x':  # y - z plane
                pos = pos[:, [1, 2, 0]]
                if self.sp is False:
                    self.cc[0], self.cc[1], self.cc[2] = s.center[1] - self.cc[1],\
                        s.center[2] - self.cc[2], s.center[0] - self.cc[0]
            else:
                if self.axis.lower() != 'z':  # project to xy plane
                    raise ValueError("Do not accept this value %s for projection" % self.axis)
                if self.sp is False:
                    self.cc = s.center - self.cc
        elif isinstance(self.axis, type([])):
            if len(self.axis) == 3:
                sa, ca = np.sin(self.axis[0] / 180. *
                                np.pi), np.cos(self.axis[0] / 180. * np.pi)
                sb, cb = np.sin(self.axis[1] / 180. *
                                np.pi), np.cos(self.axis[1] / 180. * np.pi)
                sg, cg = np.sin(self.axis[2] / 180. *
                                np.pi), np.cos(self.axis[2] / 180. * np.pi)
                # ratation matrix from
                # http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
                Rxyz = np.array(
                    [[cb * cg, cg * sa * sb - ca * sg, ca * cg * sb + sa * sg],
                     [cb * sg, ca * cg + sa * sb * sg, ca * sb * sg - cg * sa],
                     [-sb,     cb * sa,                ca * cb]], dtype=np.float64)
                pos = np.dot(pos, Rxyz)
                if self.sp is False:
                    self.cc = np.dot(s.center - self.cc, Rxyz)
            else:
                raise ValueError(
                    "Do not accept this value %s for projection" % self.axis)
        elif isinstance(self.axis, type(np.array([]))):
            if len(self.axis.shape) == 1:  # only vector from the line of sight.
                normed_axis = self.axis/np.sqrt(np.sum(self.axis**2))
                pos = np.cross(normed_axis, np.cross(pos, normed_axis))
                if self.sp is False:
                    self.cc = np.cross(normed_axis, np.cross(s.center - self.cc, normed_axis))
            elif len(self.axis.shape) == 2:
                if self.axis.shape[0] == self.axis.shape[1] == 3:
                    pos = np.dot(pos, self.axis)
                    if self.sp is False:
                        self.cc = np.dot(s.center - self.cc, self.axis)
                else:
                    raise ValueError("Axis shape is not 3x3: ", self.axis.shape)
        else:
            raise ValueError(
                "Do not accept this value %s for projection" % self.axis)

        if self.zthick is not None:
            if self.sp is False:  # pos is not centered at 0,0,0
                ids = (pos[:, 2] > self.cc[2] - self.zthick) & (pos[:, 2] < self.cc[2] + self.zthick)
            else:
                ids = (pos[:, 2] > -self.zthick) & (pos[:, 2] < self.zthick)
            pos = pos[ids]
            for i in d.keys():
                d[i] = d[i][ids]

        if self.ar is None:
            self.pxsize = np.min([pos[:, 0].max()-pos[:, 0].min(), pos[:, 1].max()-pos[:, 1].min()])/self.npx
            if self.z <= 0.0:
                self.ar = self.pxsize / s.cosmology.h * s.cosmology.arcsec_per_kpc_comoving(0.05).value
            else:
                self.ar = self.pxsize / s.cosmology.h * s.cosmology.arcsec_per_kpc_comoving(self.z).value
        else:
            if self.z <= 0.0:
                self.z = 0.05
            self.pxsize = self.ar / s.cosmology.arcsec_per_kpc_comoving(self.z).value * s.cosmology.h
        self.ar /= 3600.  # arcsec to degree

        if self.sp is not False:  # noraml projection
            minx = -(self.npx + 1) * self.pxsize / 2
            maxx = +(self.npx + 1) * self.pxsize / 2
            miny = -(self.npx + 1) * self.pxsize / 2
            maxy = +(self.npx + 1) * self.pxsize / 2
            xx = np.arange(minx, maxx, self.pxsize)
            yy = np.arange(miny, maxy, self.pxsize)
        else:  # we do real projection by transferring into RA, Dec
            SC = SkyCoord(pos[:, 0]/s.cosmology.h*u.kpc, pos[:, 1]/s.cosmology.h*u.kpc,
                          pos[:, 2]/s.cosmology.h*u.kpc, representation='cartesian')
            SC = SC.transform_to('icrs')
            pos[:, 0], pos[:, 1] = SC.ra.degree, SC.dec.degree

            SC = SkyCoord(self.cc[0]/s.cosmology.h*u.kpc, self.cc[1]/s.cosmology.h*u.kpc,
                          self.cc[2]/s.cosmology.h*u.kpc, representation='cartesian')
            SC = SC.transform_to('icrs')
            self.sp = [SC.ra.degree, SC.dec.degree]

            minx = self.sp[0] - (self.npx + 1) * self.ar / 2.
            maxx = self.sp[0] + (self.npx + 1) * self.ar / 2.
            miny = self.sp[1] - (self.npx + 1) * self.ar / 2.
            maxy = self.sp[1] + (self.npx + 1) * self.ar / 2.
            xx = np.arange(minx, maxx, self.ar)
            yy = np.arange(miny, maxy, self.ar)

        for i in d.keys():
            if self.flux:  # luminosity
                self.outd[i] = np.histogram2d(pos[:, 0], pos[:, 1], bins=[xx, yy], weights=d[i])[0]
            else:  # ab mag
                self.outd[i] = np.histogram2d(pos[:, 0], pos[:, 1], bins=[xx, yy], weights=10**(d[i]/-2.5))[0]

        # Now grid the data
        # pmax, pmin = np.max(self.S_pos, axis=0), np.min(self.S_pos, axis=0)
        # grid_x, grid_y = np.mgrid[pmin[0]:pmax[0]:nx, pmin[1]:pmax[1]:nx]
        if self.omas or self.oage or self.omet:
            self.outd["Mass"] = np.histogram2d(pos[:, 0], pos[:, 1], bins=[xx, yy], weights=s.S_mass)[0]
        if self.oage:
            ids = self.outd["Mass"] > 0
            self.outd["Age"] = np.histogram2d(pos[:, 0], pos[:, 1], bins=[xx, yy], weights=s.S_mass * s.S_age)[0]
            self.outd["Age"][ids] /= self.outd["Mass"][ids]
        if self.omet:
            ids = self.outd["Mass"] > 0
            self.outd["Metal"] = np.histogram2d(pos[:, 0], pos[:, 1], bins=[xx, yy], weights=s.S_mass * s.S_metal)[0]
            self.outd["Metal"][ids] /= self.outd["Mass"][ids]

        self.pxsize /= s.cosmology.h  # Now pixel size in physical
        self.cc = s.center/s.cosmology.h  # real center in the data

    def write_fits_image(self, fname, clobber=False):
        r"""
        Generate a image by binning X-ray counts and write it to a FITS file.

        Parameters
        ----------
        imagefile : string
            The name of the image file to write.
        clobber : boolean, optional
            Set to True to overwrite a previous file.
        """
        import pyfits as pf

        if fname[-5:] != ".fits":
            fname = fname + ".fits"

        for i in self.outd.keys():
            hdu = pf.PrimaryHDU(self.outd[i].T)
            hdu.header["SIMPLE"] = 'T'
            hdu.header.comments["SIMPLE"] = 'conforms to FITS standard'
            hdu.header["BITPIX"] = int(-32)
            hdu.header.comments["BITPIX"] = '32 bit floating point'
            hdu.header["NAXIS"] = int(2)
            hdu.header["NAXIS1"] = int(self.outd[i].shape[0])
            hdu.header["NAXIS2"] = int(self.outd[i].shape[1])
            hdu.header["EXTEND"] = True
            hdu.header.comments["EXTEND"] = 'Extensions may be present'
            hdu.header["FILTER"] = i
            hdu.header.comments["FILTER"] = 'filter used'
            hdu.header["RADECSYS"] = 'ICRS    '
            hdu.header.comments["RADECSYS"] = "International Celestial Ref. System"
            hdu.header["CTYPE1"] = 'RA---TAN'
            hdu.header.comments["CTYPE1"] = "Coordinate type"
            hdu.header["CTYPE2"] = 'DEC--TAN'
            hdu.header.comments["CTYPE2"] = "Coordinate type"
            hdu.header["CUNIT1"] = 'deg     '
            hdu.header.comments["CUNIT1"] = 'Units'
            hdu.header["CUNIT2"] = 'deg     '
            hdu.header.comments["CUNIT2"] = 'Units'
            hdu.header["CRPIX1"] = float(self.npx/2.0)
            hdu.header.comments["CRPIX1"] = 'X of reference pixel'
            hdu.header["CRPIX2"] = float(self.npx/2.0)
            hdu.header.comments["CRPIX2"] = 'Y of reference pixel'
            hdu.header["CRVAL1"] = float(self.sp[0])
            hdu.header.comments["CRVAL1"] = 'RA of reference pixel (deg)'
            hdu.header["CRVAL2"] = float(self.sp[1])
            hdu.header.comments["CRVAL2"] = 'Dec of reference pixel (deg)'
            hdu.header["CD1_1"] = float(self.ar)
            hdu.header.comments["CD1_1"] = 'RA deg per column pixel'
            hdu.header["CD1_2"] = float(0)
            hdu.header.comments["CD1_2"] = 'RA deg per row pixel'
            hdu.header["CD2_1"] = float(0)
            hdu.header.comments["CD2_1"] = 'Dec deg per column pixel'
            hdu.header["CD2_2"] = float(self.ar)
            hdu.header.comments["CD2_2"] = 'Dec deg per row pixel'

            hdu.header["RCVAL1"] = float(self.cc[0])
            hdu.header.comments["RCVAL1"] = 'Real center X of the data'
            hdu.header["RCVAL2"] = float(self.cc[1])
            hdu.header.comments["RCVAL2"] = 'Real center Y of the data'
            hdu.header["RCVAL3"] = float(self.cc[2])
            hdu.header.comments["RCVAL3"] = 'Real center Z of the data'
            hdu.header["UNITS"] = "kpc"
            hdu.header.comments["UNITS"] = 'Units for the RCVAL and PSIZE'
            hdu.header["ORAD"] = float(self.rr)
            hdu.header.comments["ORAD"] = 'Radius for cutting the object'
            hdu.header["REDSHIFT"] = float(self.z)
            hdu.header.comments["REDSHIFT"] = 'The redshift of the object'
            hdu.header["PSIZE"] = float(self.pxsize)
            hdu.header.comments["PSIZE"] = 'The pixel size of the image in comoving'
            hdu.header["AGLRES"] = float(self.ar*3600.)
            hdu.header.comments["AGLRES"] = '\'observation\' angular resolution in arcsec'
            hdu.header["ORIGIN"] = 'PymGal'
            hdu.header["VERSION"] = version.version  # get_property('__version__')
            hdu.header["DATE-OBS"] = Time.now().tt.isot
            hdu.header["NOTE"] = self.note
            hdu.writeto(fname[:-5]+"-"+i+fname[-5:], clobber=clobber)
