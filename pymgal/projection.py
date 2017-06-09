import numpy as np
# scipy must >= 0.17 to properly use this!
# from scipy.stats import binned_statistic_2d
import pyfits as pf


class projection(object):
    r"""load analysing data from simulation snapshots (gadget format only),
    yt, or raw data. Currently only works with snapshot=True. Will work more
    on other data sets.

    Parameters
    ----------
    wdata   : The data to be saved. Type: array or list of array.
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

    Notes
    -----


    Example
    -------


    """

    def __init__(self, data, simd, axis="z", npx=512, AR=None, redshift=None, flux=False):

        self.axis = axis
        self.npx = npx
        self.ar = AR
        if redshift is None:
            self.z = simd.z
        else:
            self.z = redshift
        self.flux = flux
        self.medpos = []
        self.nx = 0
        self.ny = 0

        if not flux:
            if isinstance(data, type(np.zeros(1))):
                data = 10**(data/-2.5)
                self.outd = np.zeros((npx, npx), dtype=float)
            elif isinstance(data, type([])):
                self.outd = []
                for i in len(data):
                    data[i] = 10**(data[i]/-2.5)
                    # self.outd.append(np.zeros((npx, npx), dtype=float))
            elif isinstance(data, type({})):
                self.outd = {}
                for i in data.keys():
                    data[i] = 10**(data[i]/-2.5)
                    # self.outd[i] = np.zeros((npx, npx), dtype=float)
            else:
                raise ValueError("Do not accept this data type %s " % type(data))

        if len(data) != npx**2:
            self._prep_out(self, data, simd)

    def _prep_out(self, d, s):
        r""" rotate the data points and project them into a 2D grid.
        """
        # self.nx = nx
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

        minx = pos[:, 0].min()
        maxx = pos[:, 0].max()
        miny = pos[:, 1].min()
        maxy = pos[:, 1].max()
        self.medpos = np.median(pos, axis=0)
        if self.ar is None:
            pixelsize = np.min([maxx - minx, maxy - miny]) / self.npx
        else:
            if self.z <= 0.0:
                self.z = 0.01
            pixelsize = s.cosmology.arcsec_per_kpc_proper(self.z) * self.ar * s.cosmology.h
            minx = self.medpos[0] - self.npx * pixelsize / 2
            maxx = self.medpos[0] + self.npx * pixelsize / 2
            miny = self.medpos[1] - self.npx * pixelsize / 2
            maxy = self.medpos[1] + self.npx * pixelsize / 2

        xx = np.arange(minx, maxx, pixelsize)
        self.nx = xx.size
        yy = np.arange(miny, maxy, pixelsize)
        self.ny = yy.size
        if isinstance(d, type(np.zeros(1))):
            self.outd = np.histogram2d(pos[:, 0], pos[:, 1], bins=[xx, yy], weights=d)[0]
        elif isinstance(d, type([])):
            for i in len(d):
                self.outd.append(np.histogram2d(pos[:, 0], pos[:, 1], bins=[xx, yy], weights=d[i])[0])
        elif isinstance(d, type({})):
            for i in d.keys():
                self.outd[i] = np.histogram2d(pos[:, 0], pos[:, 1], bins=[xx, yy], weights=d[i])[0]

        # Now grid the data
        # pmax, pmin = np.max(self.S_pos, axis=0), np.min(self.S_pos, axis=0)
        # grid_x, grid_y = np.mgrid[pmin[0]:pmax[0]:nx, pmin[1]:pmax[1]:nx]
        # self.grid_mass = np.histogram2d(self.S_pos[:, 0], self.S_pos[:, 1], bins=[
        #                                 nx, nx], weights=self.S_mass)[0]
        # ids = self.grid_mass > 0
        # self.grid_age = np.histogram2d(self.S_pos[:, 0], self.S_pos[:, 1], bins=[
        #                                nx, nx], weights=self.S_age * self.S_mass)[0]
        # self.grid_age[ids] /= self.grid_mass[ids]  # mass weighted age
        # self.grid_metal = np.histogram2d(self.S_pos[:, 0], self.S_pos[:, 1], bins=[
        #                                  nx, nx], weights=self.S_metal * self.S_mass)[0]
        # self.grid_metal[ids] /= self.grid_mass[ids]  # mass weighted metal
        # dx = (pmax[0] - pmin[0] + pmax[0] * 0.001) / nx
        # dy = (pmax[1] - pmin[1] + pmax[1] * 0.001) / nx
        # self.grids = np.int32(
        #     np.floor((self.S_pos[:, :2] - pmin[:2]) / np.array([dx, dy])))

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
        if fname[:-5] != ".fits":
            fname = fname + ".fits"

        if isinstance(self.outd, type(np.zeros(1))):
            hdu = pf.PrimaryHDU(self.outd.T)
            hdu.header["SIMPLE"] = "T"
            hdu.header["BITPIX"] = 8
            hdu.header["NAXIS"] = 2
            hdu.header["EXTEND"] = "T"
            hdu.header["MFORM1"] = "RA,DEC"
            hdu.header["CTYPE1"] = "RA---TAN"
            hdu.header["CTYPE2"] = "DEC--TAN"
            hdu.header["CRPIX1"] = self.nx
            hdu.header["CRPIX2"] = self.ny
            hdu.header["CRVAL1"] = float(self.medpos[0])
            hdu.header["CRVAL2"] = float(self.medpos[1])
            hdu.header["CUNIT1"] = "kpc/h"
            hdu.header["CUNIT2"] = "kpc/h"
            hdu.header["NOTE"] = ""
            hdu.writeto(fname, clobber=clobber)
        elif isinstance(self.outd, type([])):
            for i in len(self.outd):
                hdu = pf.PrimaryHDU(self.outd[i].T)
                hdu.header["SIMPLE"] = "T"
                hdu.header["BITPIX"] = 8
                hdu.header["NAXIS"] = 2
                hdu.header["EXTEND"] = "T"
                hdu.header["MFORM1"] = "RA,DEC"
                hdu.header["CTYPE1"] = "RA---TAN"
                hdu.header["CTYPE2"] = "DEC--TAN"
                hdu.header["CRPIX1"] = self.nx
                hdu.header["CRPIX2"] = self.ny
                hdu.header["CRVAL1"] = float(self.medpos[0])
                hdu.header["CRVAL2"] = float(self.medpos[1])
                hdu.header["CUNIT1"] = "kpc/h"
                hdu.header["CUNIT2"] = "kpc/h"
                hdu.header["NOTE"] = ""
                hdu.writeto(fname[:-5]+"-"+str(i)+fname[-5:], clobber=clobber)
        elif isinstance(self.outd, type({})):
            for i in self.outd.keys():
                hdu = pf.PrimaryHDU(self.outd[i].T)
                hdu.header["SIMPLE"] = "T"
                hdu.header["BITPIX"] = 8
                hdu.header["NAXIS"] = 2
                hdu.header["EXTEND"] = "T"
                hdu.header["MFORM1"] = "RA,DEC"
                hdu.header["CTYPE1"] = "RA---TAN"
                hdu.header["CTYPE2"] = "DEC--TAN"
                hdu.header["CRPIX1"] = self.nx
                hdu.header["CRPIX2"] = self.ny
                hdu.header["CRVAL1"] = float(self.medpos[0])
                hdu.header["CRVAL2"] = float(self.medpos[1])
                hdu.header["CUNIT1"] = "kpc/h"
                hdu.header["CUNIT2"] = "kpc/h"
                hdu.header["NOTE"] = ""
                hdu.writeto(fname[:-5]+"-"+i+fname[-5:], clobber=clobber)
