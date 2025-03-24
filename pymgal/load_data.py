from pymgal.readsnapsgl import readsnap
import numpy as np
from astropy.cosmology import FlatLambdaCDM
import os
from contextlib import redirect_stdout
from glob import glob
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
                 datafile=None, center=None, radius=None, read_gas=False):

        self.S_age = np.array([])
        self.S_metal = np.array([])
        self.S_mass = np.array([])
        self.S_pos = np.array([])
        self.G_pos = np.array([])  # Positions of gas particles
        self.cosmology = None  # default wmap7
        self.scale_factor = 1.0  # z = 0
        self.redshift = 0.0
        self.Uage = 0.0  # university age in Gyrs
        self.center = center
        self.radius = radius
        self.read_gas = read_gas
        self.nx = self.grid_mass = self.grid_age = self.grid_metal = None
        
        if snapshot:
            self._load_snap(snapname, nmet)
        elif yt_data is not None:
            self._load_yt(yt_data)
        elif datafile is not None:
            self._load_raw(datafile)

    def _load_snap(self, filename, nmetal):

        sample_filename = filename      # Get a sample file for testing purposes
        if os.path.isdir(filename):     # If the user provides a directory, get the first entry in that dir and set that to your sample filename
            all_files = sorted(glob(os.path.join(filename,'*')))
            sample_filename = all_files[0]  
            
        # Check if the files are in hdf5. This also needs to work if the user inputs a directory instead of a single file
        is_hdf5 = False
        if sample_filename[-4:]=='hdf5':
            is_hdf5 = True
       
        head = readsnap(filename, "HEAD", quiet=True)

        GFM_string = ""    # Some simulations (e.g. IllustrisTNG) format some quantities with "GFM_" prefix (e.g. GFM_Metallicity, GFM_SellarFormationTime)
        mass_string = "Masses"  # Many simulations write mass values as "Masses", but others (e.g. EAGLE) write it as "Mass".

        with redirect_stdout(open(os.devnull, 'w')):   # Since we are just checking for alternative names, we can suppress print statements warning that the metallicity/mass field cannot be found
            if is_hdf5 and readsnap(sample_filename, "GFM_Metallicity", ptype=4, quiet=True) is not None:   # Test if 'GFM_Metallicity' is present. If so, we will add the GFM_ prefix where needed. Otherwise leave it as is
                GFM_string = "GFM_"
            if is_hdf5 and readsnap(sample_filename, "Masses", ptype=4, quiet=True) is None and readsnap(sample_filename, "Mass", ptype=4, quiet=True) is not None:
                mass_string = "Mass"
            
        self.cosmology = FlatLambdaCDM(head.HubbleParam * 100, head.Omega0)
        self.scale_factor = head.Time
        self.redshift = head.Redshift  # if head[3] > 0 else 0.0
        self.Uage = self.cosmology.age(self.redshift).value

        
        
        if is_hdf5:  # If file is hdf5
            spos = readsnap(filename, "Coordinates", ptype=4, quiet=True)
            gpos = readsnap(filename, "Coordinates", ptype=0, quiet=True) if self.read_gas else None   # Find gas particles if requested, otherwise let g_pos be None
            self.G_mass = readsnap(filename, mass_string, ptype=0, quiet=True) if self.read_gas else None 

        else:                          # If file is not hdf5
            spos = readsnap(filename, "POS ", ptype=4, quiet=True)
            gpos = readsnap(filename, "POS ", ptype=0, quiet=True) if self.read_gas else None
            self.G_mass = readsnap(filename, "MASS", ptype=0, quiet=True) if self.read_gas else None

        
        if (self.center is not None) and (self.radius is not None):
            # r = np.sqrt(np.sum((spos - self.center)**2, axis=1))
            # ids = r <= self.radius
            ids = (spos[:, 0] >= self.center[0] - self.radius) & \
                (spos[:, 0] <= self.center[0] + self.radius) & \
                (spos[:, 1] >= self.center[1] - self.radius) & \
                (spos[:, 1] <= self.center[1] + self.radius) & \
                (spos[:, 2] >= self.center[2] - self.radius) & \
                (spos[:, 2] <= self.center[2] + self.radius)                               # Only keep stellar particles within the specified range
            self.S_pos = spos[ids]  # - self.center

            gas_ids = (gpos[:, 0] >= self.center[0] - self.radius) & \
                (gpos[:, 0] <= self.center[0] + self.radius) & \
                (gpos[:, 1] >= self.center[1] - self.radius) & \
                (gpos[:, 1] <= self.center[1] + self.radius) & \
                (gpos[:, 2] >= self.center[2] - self.radius) & \
                (gpos[:, 2] <= self.center[2] + self.radius) if self.read_gas else None    # Only keep gas particles within the specified range
            self.G_pos = gpos[gas_ids] if self.read_gas else self.G_pos
            self.center = np.asarray(self.center)
        else:
            ids = np.ones(head.totnum[4], dtype=bool)
            self.S_pos = spos  # - np.mean(spos, axis=0)
            self.G_pos = gpos if self.read_gas else self.G_pos
            self.center = np.mean(spos, axis=0)
            self.radius = np.max([spos[:,0].max()-self.center[0], self.center[0]-spos[:,0].min(),
                                  spos[:,1].max()-self.center[1], self.center[1]-spos[:,1].min(),
                                  spos[:,2].max()-self.center[2], self.center[2]-spos[:,2].min()])
        if is_hdf5:
            age = readsnap(filename, GFM_string + "StellarFormationTime",ptype=4, quiet=True)[ids]
            self.S_mass = readsnap(filename, mass_string, ptype=4, quiet=True)[ids] * 1.0e10 / head.HubbleParam   # either "Masses" or "Mass", depending on what we found above
            if head.F_Metals > 1:
                self.S_metal = readsnap(filename, GFM_string + "Metallicity", ptype=4, quiet=True)[:,0] # note this is only for SIMBA simulation
            else:
                self.S_metal = readsnap(filename, GFM_string + "Metallicity", ptype=4, quiet=True)
                
        else:
            age = readsnap(filename, "AGE ", quiet=True)[:head.totnum[4]][ids]
            self.S_mass = readsnap(filename, "MASS", ptype=4, quiet=True)[ids] * 1.0e10 / head.HubbleParam  # in M_sun
            if readsnap(filename, "Z   ", ptype=4, nmet=nmetal, quiet=True) is not None:
                self.S_metal = readsnap(filename, "Z   ", ptype=4, nmet=nmetal, quiet=True)
            else:
                self.S_metal = readsnap(filename, "ZTOT", ptype=4, nmet=nmetal, quiet=True)
            
        age = self.Uage - self.cosmology.age(1. / age - 1).value
        age[age < 0] = 0  # remove negative ages
        self.S_age = age * 1.0e9  # in yrs
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
