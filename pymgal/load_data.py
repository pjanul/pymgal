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
    nmet        : Set how many metals elements in your snapshots file, default = 11

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

    def __init__(self, snapname='', snapshot=False, nmet=11, yt_data=None,
                 datafile=None, center=None, radius=None):

        self.S_age = np.array([])
        self.S_metal = np.array([])
        self.S_mass = np.array([])
        self.S_pos = np.array([])
        self.cosmology = None  # default wmap7
        self.scale_factor = 1.0  # z = 0
        self.redshift = 0.0
        self.Uage = 0.0  # university age in Gyrs
        self.center = center
        self.radius = radius
        self.nx = self.grid_mass = self.grid_age = self.grid_metal = None

        if snapshot:
            self._load_snap(snapname, nmet)
        elif yt_data is not None:
            self._load_yt(yt_data)
        elif datafile is not None:
            self._load_raw(datafile)

    def _load_snap(self, filename, nmetal):
        head = readsnapsgl(filename, "HEAD", quiet=True)
        self.cosmology = FlatLambdaCDM(head[-1] * 100, head[-3])
        self.scale_factor = head[2]
        self.redshift = head[3] if head[3] > 0 else 0.0
        self.Uage = self.cosmology.age(self.redshift).value

        spos = readsnapsgl(filename, "POS ", ptype=4, quiet=True)
        if (self.center is not None) and (self.radius is not None):
            # r = np.sqrt(np.sum((spos - self.center)**2, axis=1))
            # ids = r <= self.radius
            ids = (spos[:, 0] >= self.center[0] - self.radius) & \
                (spos[:, 0] <= self.center[0] + self.radius) & \
                (spos[:, 1] >= self.center[1] - self.radius) & \
                (spos[:, 1] <= self.center[1] + self.radius) & \
                (spos[:, 2] >= self.center[2] - self.radius) & \
                (spos[:, 2] <= self.center[2] + self.radius)
            self.S_pos = spos[ids]  # - self.center
            self.center = np.asarray(self.center)
        else:
            ids = np.ones(head[0][4], dtype=bool)
            self.S_pos = spos  # - np.mean(spos, axis=0)
            self.center = np.mean(spos, axis=0)
            self.radius = np.max([spos[:,0].max()-self.center[0], self.center[0]-spos[:,0].min(),
                                  spos[:,1].max()-self.center[1], self.center[1]-spos[:,1].min(),
                                  spos[:,2].max()-self.center[2], self.center[2]-spos[:,2].min()])
        age = readsnapsgl(filename, "AGE ", quiet=True)[:head[0][4]][ids]
        age = self.Uage - self.cosmology.age(1. / age - 1).value
        age[age < 0] = 0  # remove negative ages
        self.S_age = age * 1.0e9  # in yrs
        self.S_mass = readsnapsgl(filename, "MASS", ptype=4, quiet=True)[
            ids] * 1.0e10 / head[-1]  # in M_sun
        self.S_metal = readsnapsgl(filename, "Z   ", ptype=4, nmet=nmetal, quiet=True)
        if self.S_metal is 0:
            self.S_metal = readsnapsgl(filename, "ZTOT", ptype=4, nmet=nmetal, quiet=True)[ids]
        else:
            self.S_metal = self.S_metal[ids]

    def _load_yt(self, yt):
        # Need to be done soon
        sp = yt.sperical(yt)
        return(sp)

    def _load_raw(self, datafile):
        if (self.center is not None) and (self.radius is not None):
            r = np.sqrt(np.sum((datafile['pos'] - self.center)**2, axis=1))
            ids = r <= self.radius
        else:
            ids = np.ones(datafile['age'].size, dtype=bool)
        self.S_age = datafile['age'][ids]
        self.S_mass = datafile['mass'][ids]
        self.S_metal = datafile['metal'][ids]
