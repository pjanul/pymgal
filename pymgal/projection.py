import numpy as np
# scipy must >= 0.17 to properly use this!
# from scipy.stats import binned_statistic_2d


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
    axis    : can be 'x', 'y', 'z', or a list of degrees [alpha, beta, gamma],
                which will rotate the data points by $\alpha$ around the x-axis,
                $\beta$ around the y-axis, and $\gamma$ around the z-axis
                Default: "z"
    npx     : The pixel number of the grid for projection.
                Type: int. Default: 512
                A npx x npx image will be produced later.
    AR      : Angular resolution. Type: arcsec. Default: None
    redshift: The redshift of the object at. Default: None.
                If None, redshift from simulation data will be used.
                Maybe move to 0.01 if it is 0.
    flux    : is the input wdata in luminosity/flux? Default: False, assume in ab mag.
                Set this to true if particles' luminosity are given in wdata.
    outmas  : do you want to out put stellar mass? Default: False.
                If True, the stellar mass in each pixel are saved.
    outage  : do you want to out put stellar age (mass weighted)? Default: False.
                If True, the stellar age in each pixel are saved.
    outmet  : do you want to out put stellar metallicity (mass weighted)? Default: False.
                If True, the stellar metallicity in each pixel are saved.
    Notes
    -----


    Example
    -------
    Pdata = pymgal.projection(part_lum, simu_data, npx=1024, flux=Ture)
    Pdata.write_fits_image("filename.fits")
    """

    def __init__(self, data, simd, axis="z", npx=512, AR=None, redshift=None,
                 flux=False, outmas=False, outage=False, outmet=False):

        self.axis = axis
        self.npx = npx
        self.ar = AR
        if redshift is None:
            self.z = simd.z
        else:
            self.z = redshift
        self.flux = flux
        self.nx = 0
        self.ny = 0
        self.pxsize = 0.
        self.cc = simd.center
        self.rr = simd.radius
        self.flux = flux
        self.omas = outmas
        self.oage = outage
        self.omet = outmet
        self.outd = {}

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

        if isinstance(self.axis, type('')):
            pos = s.S_pos[:, :2]
            if self.axis.lower() == 'y':  # x-z plane
                pos[:, 1] = s.S_pos[:, 2]
            elif self.axis.lower() == 'x':  # y - z plane
                pos = s.S_pos[:, 1:]
            else:
                if self.axis.lower() != 'z':  # project to xy plane
                    raise ValueError(
                        "Do not accept this value %s for projection" % self.axis)
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
                pos = np.dot(s.S_pos, Rxyz)[:, :2]
            else:
                raise ValueError(
                    "Do not accept this value %s for projection" % self.axis)
        else:
            raise ValueError(
                "Do not accept this value %s for projection" % self.axis)

        minx = pos[:, 0].min()
        maxx = pos[:, 0].max()
        miny = pos[:, 1].min()
        maxy = pos[:, 1].max()
        if self.ar is None:
            self.pxsize = np.min([maxx - minx, maxy - miny]) / self.npx
        else:
            if self.z <= 0.0:
                self.z = 0.01
            self.pxsize = self.ar / s.cosmology.arcsec_per_kpc_proper(self.z).value * s.cosmology.h
            minx = -self.npx * self.pxsize / 2
            maxx = +self.npx * self.pxsize / 2
            miny = -self.npx * self.pxsize / 2
            maxy = +self.npx * self.pxsize / 2

        xx = np.arange(minx, maxx, self.pxsize)
        self.nx = xx.size
        yy = np.arange(miny, maxy, self.pxsize)
        self.ny = yy.size

        for i in d.keys():
            if self.flux:  # luminosity
                self.outd[i] = np.histogram2d(pos[:, 0], pos[:, 1], bins=[xx, yy], weights=d[i])[0]
            else:  # ab mag
                self.outd[i] = np.histogram2d(pos[:, 0], pos[:, 1],
                                              bins=[xx, yy], weights=10**(d[i]/-2.5))[0]

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
            hdu.header["RCVAL1"] = float(self.cc[0])
            hdu.header["RCVAL2"] = float(self.cc[1])
            hdu.header["UNITS"] = "kpc/h"
            hdu.header["ORAD"] = float(self.rr)
            hdu.header["REDSHIFT"] = float(self.z)
            hdu.header["PSIZE"] = float(self.pxsize)
            hdu.header["AGLRES"] = float(self.ar)
            hdu.header["NOTE"] = ""
            hdu.writeto(fname[:-5]+"-"+i+fname[-5:], clobber=clobber)
