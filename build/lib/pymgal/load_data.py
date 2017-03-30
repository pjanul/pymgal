from pymgal.readsnapsgl import readsnapsgl
import numpy as np
from astropy.cosmology import FlatLambdaCDM
# from scipy.interpolate import griddata


class load_data(object):
    r"""load analysing data from simulation snapshots (gadget format only),
    yt, or raw data. Currently only works with snapshot=True. Will work more
    on other data sets.

    Parameters
    ----------
    snapname    : The filename of simulation snapshot. Default : ''
    snapshot    : Is loading snapshot or not? Default : False

    yt_data     : yt data set ds (ds = yt.load()). Default : None
                  Don't use! Currently do not support yet.

    datafile    : raw data set, can be file, np array, or dictionary. Default : None
                  Don't use! Currently do not support yet.

    ----------If you only need parts of loaded data. Specify center and radius
    center      : The center of a sphere for the data you want to get.
                  Default : None
    radius      : The radius of a sphere for the data you want to get.
                  Default : None

    Notes
    -----
    Need to take more care about the units!!! Currently only assume simulation units.
    kpc/h and 10^10 M_sun
    Raw data set needs to provide the cosmology ...

    Example
    -------
    simd = load_data(snapfilename="/home/weiguang/Downloads/snap_127",
                     snapshot=True,
                     center=[500000,500000,500000], radius=800)

    """

    def __init__(self, snapname='', snapshot=False, yt_data=None,
                 datafile=None, center=None, radius=None):

        self.S_age = np.array([])
        self.S_metal = np.array([])
        self.S_mass = np.array([])
        self.S_pos = np.array([])
        self.cosmology = None  # default wmap7
        self.currenta = 1.0  # z = 0
        self.z = 0.0
        self.Uage = 0.0  # university age in Gyrs
        self.nx = self.grid_mass = self.grid_age = self.grid_metal = None

        if snapshot:
            self._load_snap(snapname, center, radius)
        elif yt_data is not None:
            self._load_yt(yt_data, center, radius)
        elif datafile is not None:
            self._load_raw(datafile, center, radius)

    def _load_snap(self, filename, cc, rr):
        head = readsnapsgl(filename, "HEAD", quiet=True)
        self.cosmology = FlatLambdaCDM(head[-1] * 100, head[-3])
        self.currenta = head[2]
        self.Uage = self.cosmology.age(1. / self.currenta - 1)
        self.z = head[3] if head[3] > 0 else 0.0
        spos = readsnapsgl(filename, "POS ", ptype=4, quiet=True)
        if (cc is not None) and (rr is not None):
            r = np.sqrt(np.sum((spos - cc)**2, axis=1))
            ids = r <= rr
            self.S_pos = spos[ids] - cc
        else:
            ids = np.ones(head[0][4], dtype=bool)
            self.S_pos = spos - np.mean(spos, axis=0)
        age = readsnapsgl(filename, "AGE ", ptype=4, quiet=True)[ids]
        age = self.Uage - self.cosmology.age(1. / age - 1)
        self.S_age = age.value * 1.0e9  # in yrs
        self.S_mass = readsnapsgl(filename, "MASS", ptype=4, quiet=True)[
            ids] * 1.0e10 / head[-1]  # in M_sun
        self.S_metal = readsnapsgl(filename, "Z   ", ptype=4, quiet=True)[ids]

    def _load_yt(self, yt, cc, rr):
        # Need to be done soon
        sp = yt.sperical(yt)
        return(sp)

    def _load_raw(self, datafile, cc, rr):
        if (cc is not None) and (rr is not None):
            r = np.sqrt(np.sum((datafile['pos'] - cc)**2, axis=1))
            ids = r <= rr
        else:
            ids = np.ones(datafile['age'].size, dtype=bool)
        self.S_age = datafile['age'][ids]
        self.S_mass = datafile['mass'][ids]
        self.S_metal = datafile['metal'][ids]

    def rotate_grid(self, axis, nx):
        r""" rotate the data points and project them into a 2D grid.

        Parameter:
        ----------
        axis    : can be 'x', 'y', 'z', or a list of degrees [alpha, beta, gamma],
                  which will rotate the data points by $\alpha$ around the x-axis,
                  $\beta$ around the y-axis, and $\gamma$ around the z-axis
        nx      : The pixel size of the grid. A nx x nx image can be produced later.
        """
        self.nx = nx
        # ratation data points first
        if isinstance(axis, type('')):
            if axis == 'y':  # x-z plane
                self.S_pos[:, 1] = self.S_pos[:, 2]
            elif axis == 'x':  # y - z plane
                self.S_pos[:, 0] = self.S_pos[:, 1]
                self.S_pos[:, 1] = self.S_pos[:, 2]
            else:
                if axis != 'z':  # project to xy plane
                    raise ValueError(
                        "Do not accept this value %s for projection" % axis)
        elif isinstance(axis, type([])):
            if len(axis) == 3:
                sa, ca = np.sin(axis[0] / 180. *
                                np.pi), np.cos(axis[0] / 180. * np.pi)
                sb, cb = np.sin(axis[1] / 180. *
                                np.pi), np.cos(axis[1] / 180. * np.pi)
                sg, cg = np.sin(axis[2] / 180. *
                                np.pi), np.cos(axis[2] / 180. * np.pi)
                # ratation matrix from
                # http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
                Rxyz = np.array(
                    [[cb * cg, cg * sa * sb - ca * sg, ca * cg * sb + sa * sg],
                     [cb * sg, ca * cg + sa * sb * sg, ca * sb * sg - cg * sa],
                     [-sb,     cb * sa,                ca * cb]], dtype=np.float64)
                self.S_pos = np.dot(self.S_pos, Rxyz)
            else:
                raise ValueError(
                    "Do not accept this value %s for projection" % axis)

        # Now grid the data
        # pmax, pmin = np.max(self.S_pos, axis=0), np.min(self.S_pos, axis=0)
        # grid_x, grid_y = np.mgrid[pmin[0]:pmax[0]:nx, pmin[1]:pmax[1]:nx]
        self.grid_mass = np.histogram2d(self.S_pos[:, 0], self.S_pos[:, 1], bins=[
                                        nx, nx], weights=self.S_mass)[0]
        ids = self.grid_mass > 0
        self.grid_age = np.histogram2d(self.S_pos[:, 0], self.S_pos[:, 1], bins=[
                                       nx, nx], weights=self.S_age * self.S_mass)[0]
        self.grid_age[ids] /= self.grid_mass[ids]  # mass weighted age
        self.grid_metal = np.histogram2d(self.S_pos[:, 0], self.S_pos[:, 1], bins=[
                                         nx, nx], weights=self.S_metal * self.S_mass)[0]
        self.grid_metal[ids] /= self.grid_mass[ids]  # mass weighted metal
        # dx = (pmax[0] - pmin[0] + pmax[0] * 0.001) / nx
        # dy = (pmax[1] - pmin[1] + pmax[1] * 0.001) / nx
        # self.grids = np.int32(
        #     np.floor((self.S_pos[:, :2] - pmin[:2]) / np.array([dx, dy])))
