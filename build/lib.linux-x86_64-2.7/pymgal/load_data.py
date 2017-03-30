from pymgal.readsnapsgl import readsnapsgl
import numpy as np
from astropy.cosmology import FlatLambdaCDM


class load_data(object):
    r"""
    load_data can load simulation data from snapshots (gadget format), yt, or raw data.
    Currently only works with snapshot=True
    """

    def __init__(self, snapfilename='', snapshot=False, yt_data=None,
                 datafile=None, center=None, radius=None):

        self.S_age = np.array([])
        self.S_metal = np.array([])
        self.S_mass = np.array([])
        self.S_pos = np.array([])
        self.cosmology = None  # default wmap7
        self.currenta = 1.0  # z = 0
        self.z = 0.0
        self.Uage = 0.0  # university age in Gyrs

        if snapshot:
            self._load_snap(snapfilename, center, radius)
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
        pmax, pmin = np.max(self.S_pos, axis=0), np.min(self.S_pos, axis=0)
        grid_x, grid_y = np.mgrid[pmin[0]:pmax[0]:nx, pmin[1]:pmax[1]:nx]
        self.grid_mass = griddata(
            self.S_pos[:, :2], self.S_mass, (grid_x, grid_y))
        self.grid_age = griddata(self.S_pos[:, :2], self.S_age * self.S_mass,
                                 (grid_x, grid_y)) / grid_mass  # mass weighted age
        self.grid_metal = griddata(self.S_pos[:, :2], self.S_metal * self.mass,
                                   (grid_x, grid_y)) / grid_mass  # mass weighted metal
        dx = (pmax[0] - pmin[0] + pmax[0] * 0.001) / nx
        dy = (pmax[1] - pmin[1] + pmax[1] * 0.001) / nx
        self.grids = np.int32(
            np.floor((self.S_pos[:, :2] - pmin[:2]) / np.array([dx, dy])))
