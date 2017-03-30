from readsnapsgl import readsnapsgl
import numpy as np
from astropy.cosmology import FlatLambdaCDM


class load_data(object):
    r"""
    load_data can load simulation data from snapshots (gadget format), yt, or raw data.
    Currently only works with snapshot=True
    """

    def __init__(self, snapfilename='', snapshot=False, yt_data=None, datafile=None, center=None, radius=None):

        self.S_age = np.array([])
        self.S_metal = np.array([])
        self.S_mass = np.array([])
        self.cosmology = None  # default wmap7
        self.currenta = 1.0  # z = 0
        self.z = 0.0

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
        self.z = head[3] if head[3] > 0 else 0.0
        spos = readsnapsgl(filename, "POS ", ptype=4, quiet=True)
        if (cc is not None) and (rr is not None):
            r = np.sqrt(np.sum((spos - cc)**2, axis=1))
            ids = r <= rr
        else:
            ids = np.ones(head[0][4], dtype=bool)
        age = readsnapsgl(filename, "AGE ", ptype=4, quiet=True)[ids]
        age = self.cosmology.age(1. / self.currenta) - \
            self.cosmology.age(1. / age)
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
