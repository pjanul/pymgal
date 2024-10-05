import numpy as np
import astropy.units as u
import re
from astropy.coordinates import SkyCoord, ICRS
from astropy.time import Time
from functools import lru_cache
from pymgal import __version__
# scipy must >= 0.17 to properly use this!
# from scipy.stats import binned_statistic_2d
from numba import njit, prange
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pymgal import utils



# Cache the kernels to avoid recomputation and significantly improve runtime
@lru_cache(maxsize=256)
def create_gaussian_kernel(sigma):
    x = np.arange(-3 * sigma, 3 * sigma + 1, 1)
    y = np.arange(-3 * sigma, 3 * sigma + 1, 1)
    # xx_kernel, yy_kernel = np.meshgrid(x, y)
    # kernel = np.exp(-(xx_kernel**2 + yy_kernel**2) / (2.0 * sigma**2)) / (2.0 * np.pi * sigma**2)
    xx = x[:, np.newaxis]  # Reshape x to be a column vector
    yy = y[np.newaxis, :]  # Reshape y to be a row vector
    
    kernel = np.exp(-(xx**2 + yy**2) / (2.0 * sigma**2))  / (2.0 * np.pi * sigma**2)
   
    return kernel


def add_array_to_large(large_array, small_array, center_x, center_y):
    small_channels, small_height, small_width = small_array.shape
    large_channels, large_height, large_width = large_array.shape

    # Define the boundaries of the large array and then the small array
    start_x = max(0, center_x - small_width // 2)
    end_x = min(large_width, center_x + (small_width + 1) // 2)
    start_y = max(0, center_y - small_height // 2)
    end_y = min(large_height, center_y + (small_height + 1) // 2)

    small_start_x = max(0, small_width // 2 - center_x)
    small_end_x = small_start_x + (end_x - start_x)
    small_start_y = max(0, small_height // 2 - center_y)
    small_end_y = small_start_y + (end_y - start_y)

    # Add the small array to the large one
    large_array[:, start_x:end_x, start_y:end_y] += small_array[:, small_start_x:small_end_x, small_start_y:small_end_y]
    return large_array

def gaussian_smoothing(distances, weights_dict, sample_hist, x_bins, y_bins, pixel_scale, max_kernel=None):
    
    num_channels = len(weights_dict)
    smoothed_hists = np.zeros((num_channels,) + sample_hist.shape)
    x_bins_range, y_bins_range = max(x_bins) - min(x_bins), max(y_bins) - min(y_bins)
    
    # Remove indices on image boundaries since this causes strange edge effects
    valid_indices = np.logical_and.reduce((x_bins > 0.01 * x_bins_range, x_bins < 0.99 * x_bins_range, y_bins > 0.01 * y_bins_range,y_bins < 0.99 * y_bins_range))
    distances = np.array(distances)[valid_indices]
    weights_dict = {channel: weights[valid_indices] for channel, weights in weights_dict.items()}
    x_bins, y_bins = x_bins[valid_indices], y_bins[valid_indices]
    
    # Convert the dictionary into a multi-channel array
    weights_array = np.column_stack([weights for weights in weights_dict.values()])
    
    # Calculate max sigma and make sure it is at least 1 and at most 10% of the total image dimension to keep runtime manageable 
    if max_kernel is None:
        max_kernel = int(np.round(x_bins_range / 10))
    sigmas = np.round(distances / pixel_scale)
    sigmas = np.clip(sigmas, 1, max_kernel)  # Cap sigmas between 1 and max_kernel

    # Compute the n-channel kernel for each particle and add it to the final result
    for i in range(len(distances)):
        kernel = create_gaussian_kernel(sigmas[i])
        kernel = np.stack([kernel] * num_channels, axis=0)
        kernel *= weights_array[i][:, np.newaxis, np.newaxis]
        smoothed_hists = add_array_to_large(smoothed_hists, kernel, x_bins[i], y_bins[i])

    smoothed_hists = {channel: smoothed_hists[i] for i, (channel, weights) in enumerate(weights_dict.items())}
    return smoothed_hists



def get_property(prop):
    result = re.search(r'{}\s*=\s*[\'"]([^\'"]*)[\'"]'.format(prop), open('__init__.py').read())
    return result.group(1)


class projection(object):
    r"""load analysing data from simulation snapshots,
    yt, or raw data. Currently only works with snapshot=True. May be adapted to other datasets in the future.

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
                Default: "z"
    npx     : The pixel number of the grid for projection.
                Type: int. Default: 512
                A [npx, npx] image will be produced later.
                It accept 'auto' parameter, which will automatically decide the pixel number to include all particles.
    AR      : Angular resolution. Type: arcsec. Default: None
                If AR is None and npx is auto, we will force the npx to be 512
                If AR is None, it will automatically recalculated in the code for output.
                At z=0, no matter AR is set or not, AR is always assume at z = 0.05.
    redshift: The redshift of the object at. Default: None.
                If None, redshift from simulation data will be used.
                This will be moved to 0.05 if simulation redshift is used and equal to 0.
                Note this redshift does not change anything of the simulation physical particle positions, only shift the object to this redshift for observing.
    p_thick  : The thickness in projection direction. Default: None.
                If None, use all data from cutting region. Otherwise set a value in simulation
                length unit (kpc/h normally), then a slice of data [center-p_thick, center+p_thick]
                will be used to make the y-map.
    SP      : Faked sky positions in [RA (longitude), DEC (latitude)] in degrees.
                Default: None, calculate automatically based on 3D position.  
                If [x,y,z] (len(SP) == 3) of the Earth position in the pos coordinate is given,  
                The pos - [x,y,z] are taken as the J2000 3D coordinates and converted into RA, DEC.  
    unit    : is the input wdata in luminosity, flux, Fv, Fl, Jy, magnitude? Default: flux.
                Set this to true if particles' luminosity are given in wdata.
    mag_type: If the user selected unit="magnitude", keep track of the type of magnitude (AB, vega, solar) in either apparent or absolute. Default: empty string "".
    ksmooth : An integer representing the k in kNN Gaussian smoothing. 1 sigma for the Gaussian is set to the distance to the kth neighbour in 3D space.
                Recommended value: somewhere between 20 and 80.
                If k>0, you set the smoothing length of a particle to be the distance between the particle and its kth nearest neighbour.
                If k=0, pymgal does not perform smoothing.
                If k<0, throw an error.
    lsmooth : An array of floats where each float is the smoothing length (ie kNN distance) for a given particle.
                Default: None, but can be set to a precomputed array to avoid redundant calculations.
    g_soft  : The gravitational softening length in kpc (physical), which is coverted into pixel values and used as the maximum number of pixels used for the gaussian kernel's 1 sigma.
    noise   : The noise level in Lsun/arcsec^2
    outmas  : do you want to output stellar mass? Default: False.
                If True, the stellar mass in each pixel are saved.
    outage  : do you want to out put stellar age (mass weighted)? Default: False.
                If True, the stellar age in each pixel are saved.
    outmet  : do you want to out put stellar metallicity (mass weighted)? Default: False.
                If True, the stellar metallicity in each pixel are saved.
    Notes
    -----


    Example
    -------
    Pdata = pymgal.projection(part_lum, simu_data, npx=1024, unit='flux')
    Pdata.write_fits_image("filename.fits")
    """

    def __init__(self, data, simd, axis="z", npx=512, AR=None, redshift=None, p_thick=None,
                 SP=None, unit='flux', mag_type="", ksmooth=0, lsmooth=None, g_soft=5, noise=0, outmas=False, outage=False, outmet=False):

        self.axis = axis
        if isinstance(npx, type("")) or isinstance(npx, type('')):
            self.npx = npx.lower()
        else:
            self.npx = npx
        self.ar = AR
        if redshift is None:
            self.z = simd.redshift
        else:
            self.z = redshift
        self.pxsize = 0.
        self.sp = SP
        self.cc = simd.center    
        self.rr = simd.radius / simd.cosmology.h / (1.+ simd.redshift)   # to physical in simulation time
        self.flux = unit
        self.mag_type = mag_type
        self.g_soft = g_soft / simd.cosmology.h / (1.+ simd.redshift) if g_soft is not None else None  # to physical in simulation time
        self.omas = outmas
        self.oage = outage
        self.omet = outmet
        self.ksmooth = ksmooth
        self.noise = noise
        if ksmooth < 0:
            raise ValueError("ksmooth should be a non-negative integer")
        if g_soft is not None and g_soft <= 0:
            raise ValueError("g_soft should be strictly greater than zero")
        self.lsmooth = lsmooth
        self.p_thick = p_thick
        if p_thick is not None:
            self.p_thick /= (simd.cosmology.h / (1.+ simd.redshift))
        self.outd = {}

        self._prep_out(data, simd)

    def _prep_out(self, d, s):
        r""" rotate the data points and project them into a 2D grid.
        """
            
        pos = np.copy(s.S_pos) / s.cosmology.h / (1.+ s.redshift)  # to assumed physical
        pos -= self.cc / s.cosmology.h / (1.+ s.redshift)
        center = s.center  / s.cosmology.h / (1.+ s.redshift)
        kth_distance = 0  # initialize the variable and update it if self.ksmooth > 0
        dt = d.copy()    # Save a copy of the dict and use that so you don't mess with the original one
        masses = s.S_mass
        ages = s.S_age
        metals = s.S_metal

        # If you're given a projection vector array, convert it to an angle
        if isinstance(self.axis, type(np.array([]))) and len(self.axis.shape) == 1 and len(self.axis) == 3:
            vector = np.array(self.axis, dtype=float)
            vector /= np.linalg.norm(vector)
            
            angles_list = [np.degrees(np.arctan2(vector[1], vector[2])), np.degrees(np.arctan2(-vector[0]**2,  vector[1]**2 + vector[2]**2)), 0]
            self.axis = angles_list
        
        if isinstance(self.axis, type('')):
            if self.axis.lower() == 'y':  # x-z plane
                pos = pos[:, [0, 2, 1]]
            elif self.axis.lower() == 'x':  # y - z plane
                pos = pos[:, [2, 1, 0]]
            else:
                if self.axis.lower() != 'z':  # project to xy plane
                    raise ValueError("Do not accept this value %s for projection" % self.axis)
        elif isinstance(self.axis, type([])):
            if len(self.axis) == 3:
                sa, ca = np.sin(self.axis[0] / 180. *
                                np.pi), np.cos(self.axis[0] / 180. * np.pi)
                sb, cb = np.sin(self.axis[1] / 180. *
                                np.pi), np.cos(self.axis[1] / 180. * np.pi)
                sg, cg = np.sin(self.axis[2] / 180. *
                                np.pi), np.cos(self.axis[2] / 180. * np.pi)
                # rotation matrix from
                # http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
                Rxyz = np.array(
                    [[cb * cg, cg * sa * sb - ca * sg, ca * cg * sb + sa * sg],
                     [cb * sg, ca * cg + sa * sb * sg, ca * sb * sg - cg * sa],
                     [-sb,     cb * sa,                ca * cb]], dtype=np.float64)
                pos = np.dot(pos, Rxyz)
            else:
                raise ValueError(
                    "Do not accept this value %s for projection" % self.axis)
        else:
            raise ValueError(
                "Do not accept this value %s for projection" % self.axis)
        lsmooth = self.lsmooth
        if self.p_thick is not None:
            ids = (pos[:, 2] > -self.p_thick) & (pos[:, 2] < self.p_thick)
            old_pos_size = len(pos) # keep track of the number of particles we had before we cut the z thick 
            pos = pos[ids]
            lsmooth = lsmooth[ids] if lsmooth is not None else None
            masses = masses[ids]
            ages = ages[ids]
            metals = metals[ids]
            for i, vals in d.items():
                dt[i] = d[i][ids]
        if self.ar is None:
            if self.npx == 'auto':
                self.npx = 512
            self.pxsize = np.min([pos[:, 0].max()-pos[:, 0].min(), pos[:, 1].max()-pos[:, 1].min()])/self.npx

            if self.z <= 0.0:
                self.ar = self.pxsize * s.cosmology.arcsec_per_kpc_proper(0.05).value
            else:
                self.ar = self.pxsize * s.cosmology.arcsec_per_kpc_proper(self.z).value
        else:
            if self.z <= 0.0:
                self.z = 0.05
            self.pxsize = self.ar / s.cosmology.arcsec_per_kpc_proper(self.z).value
            if self.npx == 'auto':
                self.npx = np.int32(2. * self.rr / self.pxsize) + 1
            else:
                self.rr = (self.npx-1) * self.pxsize / 2.
        self.ar /= 3600.  # arcsec to degree

        minx = -(self.npx + 1) * self.pxsize / 2
        maxx = +(self.npx + 1) * self.pxsize / 2
        miny = -(self.npx + 1) * self.pxsize / 2
        maxy = +(self.npx + 1) * self.pxsize / 2
        xx = np.arange(minx, maxx, self.pxsize)
        yy = np.arange(miny, maxy, self.pxsize)

        SC = SkyCoord(self.cc[0] * u.kpc, self.cc[1] * u.kpc, self.cc[2] * u.kpc,
                          representation_type='cartesian') 
        SC = SC.spherical
        self.sp = [SC.lon.degree, SC.lat.degree]   # longitude and latitude are equivalent to RA, dec

        x_bins = np.digitize(pos[:, 0], xx)                                                     # define x and y bins
        y_bins = np.digitize(pos[:, 1], yy)


        
        #self.cc = center  # real center in the data
        # If smoothing is set to off
        if self.ksmooth == 0:
            for i in dt.keys():
                # If magnitude, we need to remove spots with zero flux and replace them with small values
                if self.flux.lower() == "magnitude":
                    hist = np.histogram2d(pos[:, 0], pos[:, 1], bins=[xx, yy], weights=10**(-dt[i]/2.5))[0]
                    min_non_zero = np.min(hist[hist > 0])
                    hist[hist == 0] = min_non_zero / 2.0
                    self.outd[i] = -2.5*np.log10(hist)

                else: 
                    self.outd[i] = np.histogram2d(pos[:, 0], pos[:, 1], bins=[xx, yy], weights=dt[i])[0]

            if self.omas or self.oage or self.omet:
                mass_hist = np.histogram2d(pos[:, 0], pos[:, 1], bins=[xx, yy], weights=masses)[0]
                ids = mass_hist > 0

                if self.omas:
                    self.outd["Mass"] = mass_hist
                if self.oage:
                    self.outd["Age"] = np.histogram2d(pos[:, 0], pos[:, 1], bins=[xx, yy], weights=masses * ages)[0]
                    self.outd["Age"][ids] /= mass_hist[ids]
                if self.omet:
                    self.outd["Metal"] = np.histogram2d(pos[:, 0], pos[:, 1], bins=[xx, yy], weights=masses * metals)[0]
                    self.outd["Metal"][ids] /= mass_hist[ids]
                    
        # If smoothing is set to on and we have at least one filter to smooth
        elif (self.ksmooth > 0):  #and (len(d.keys()) > 0): #user is already forced to include at least one filter
            sample_hist =  np.histogram2d(pos[:, 0], pos[:, 1], bins=[xx, yy], weights=list(dt.values())[0])[0] # define a sample hist using the weights from whichever filter comes first
            smoothed_hists = None # initialize the set of smoothed histograms 
            if self.flux.lower() == "magnitude": # ab mag
                d_val = {key: 10**(-dt[key]/2.5) for key in dt}
                smoothed_hists = gaussian_smoothing(lsmooth, d_val, sample_hist, x_bins, y_bins, self.pxsize)
            else:
                d_val = {key: dt[key] for key in dt}
                smoothed_hists = gaussian_smoothing(lsmooth, d_val, sample_hist, x_bins, y_bins, self.pxsize)
            # Set the histograms to match the filters
            for i in dt.keys():
                if self.flux.lower() == 'magnitude':
                    hist = smoothed_hists[i]
                    min_non_zero = np.min(hist[hist > 0])
                    hist[hist == 0] = min_non_zero / 2.0
                    self.outd[i] = -2.5*np.log10(hist)
                else:
                    self.outd[i] = smoothed_hists[i]


            if self.omas or self.oage or self.omet:
                max_sigma = self.g_soft/self.pxsize if self.g_soft is not None else None # Set the max standard deviation for the smoothing gaussian to be no more than the gravitational softening length of 5 kpc/h
                unsmooth_mass_hist = np.histogram2d(pos[:, 0], pos[:, 1], bins=[xx, yy], weights=masses)[0]
                #ids = mass_hist >0
                hist_dict = gaussian_smoothing(lsmooth, {"Mass": masses, "Mass_Age": masses * ages, "Mass_Metal": masses * metals}, unsmooth_mass_hist, x_bins, y_bins, self.pxsize, max_kernel=max_sigma)
                mass_hist = hist_dict["Mass"]
                ids = mass_hist > 0
                if self.omas:
                    self.outd["Mass"] = mass_hist
                if self.oage:
                    age_hist = hist_dict["Mass_Age"]
                    age_hist[ids] /= mass_hist[ids]
                    self.outd["Age"] = age_hist
                if self.omet:
                    met_hist = hist_dict["Mass_Metal"]
                    met_hist[ids] /= mass_hist[ids]  
                    self.outd["Metal"] =  met_hist      

        pixel_area = (self.ar * 3600) ** 2 # Convert AR back to units from deg and compute the area of a pixel in arcsec^2
        #print(pixel_area, self.noise * pixel_area)
        for i in dt.keys():
            noise_stdev = self.noise[i]
            if noise_stdev is None:
                continue
            if self.flux == "magnitude":
                print("self noise sdtev: ", noise_stdev) 
                lum_noise = 10**(noise_stdev/-2.5)
                print(lum_noise)
                lum_vals = 10**(self.outd[i]/-2.5)
                
                lum_vals += np.random.normal(loc=0.0, scale=lum_noise * pixel_area, size=self.outd[i].shape)
                #min_non_zero = np.min(lum_vals[lum_vals > 0])
                #lum_vals[lum_vals <= 0] = min_non_zero / np.abs(np.random.normal(loc=0.0, scale=2)) # randomly add a little bit to the minimum values

                self.outd[i] = -2.5 * np.log10(lum_vals) 
            else:
                print("self.noise[i]", i, noise_stdev)
                self.outd[i] += np.random.normal(loc=0.0, scale=noise_stdev * pixel_area, size=self.outd[i].shape)

    def write_fits_image(self, fname, comments='None', overwrite=False):
        r"""
        Generate a image by binning X-ray counts and write it to a FITS file.

        Parameters
        ----------
        imagefile : string
            The name of the image file to write.
        overwrite : boolean, optional
            Set to True to overwrite a previous file.
        comments  : The comments in str will be put into the fit file header. Defualt: 'None'
                    It accepts str or list of str or tuple of str
        """
        import astropy.io.fits as pf

        if fname[-5:] != ".fits":
            fname = fname + ".fits"

        for i in self.outd.keys():
            hdu = pf.PrimaryHDU(self.outd[i].T)
            hdu.header["SIMPLE"] = True
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
            hdu.header["CD1_1"] = -float(self.ar)
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
            hdu.header["PIXVAL"] = "years" if i.lower() == "age" else "metallicity" if i.lower() == "metal" else "M_sun" if i.lower() == "mass" \
                                    else "erg/s" if self.flux == "luminosity" else "erg/s/cm^2" if self.flux == "flux" else "erg/s/cm^2/Hz" if self.flux == "fv" \
                                    else "jansky" if self.flux == "jy" else "erg/s/cm^2/angstrom" if self.flux == "fl" else "mag" + " (" + self.mag_type + ")"
            hdu.header.comments["PIXVAL"] = 'The units of the pixel values in your image'
            hdu.header["ORAD"] = float(self.rr)
            hdu.header.comments["ORAD"] = 'Rcut in physical for the image.'
            hdu.header["REDSHIFT"] = float(self.z)
            hdu.header.comments["REDSHIFT"] = 'The redshift of the object being put to'
            hdu.header["PSIZE"] = float(self.pxsize)
            hdu.header.comments["PSIZE"] = 'The pixel size in physical at simulation time'
            hdu.header["AGLRES"] = float(self.ar*3600.)
            hdu.header.comments["AGLRES"] = '\'observation\' angular resolution in arcsec'
            hdu.header["P_THICK"] = float(self.p_thick) if self.p_thick is not None else 'None'
            hdu.header.comments["P_THICK"] = 'The thickness of the projection (kpc)'
            hdu.header["ORIGIN"] = 'PyMGal'
            hdu.header.comments["ORIGIN"] = 'Software for generating this mock image'
            hdu.header["VERSION"] = __version__.__version__
            hdu.header.comments["VERSION"] = 'Version of the software'
            hdu.header["DATE-OBS"] = Time.now().tt.isot
            if isinstance(comments, type([])) or isinstance(comments, type(())):
                for j in range(len(comments)):
                    hdu.header["COMMENT"+str(j+1)] = comments[j]
            elif isinstance(comments, type("")) or isinstance(comments, type('')):
                hdu.header["COMMENT"] = comments
            else:
                raise ValueError("Do not accept this comments type! Please use str or list")
            hdu.writeto(fname[:-5]+"-"+i+fname[-5:], overwrite=overwrite)