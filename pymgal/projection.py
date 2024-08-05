import numpy as np
import astropy.units as u
import re
from astropy.coordinates import SkyCoord
from astropy.time import Time
from scipy.spatial import KDTree
import pandas as pd
from functools import lru_cache
from pymgal import __version__
# scipy must >= 0.17 to properly use this!
# from scipy.stats import binned_statistic_2d


# Create kD tree and extract the distance of the kth nearest neighbour (kth distance = 1 sigma of gaussian kernel)
def knn_distance(positions, k):
    kdtree = KDTree(positions) 
    _, indices = kdtree.query(positions, k + 1)  # k + 1 to exclude the point itself
    kth_neighbors = np.take(positions, indices[:, k], axis=0)
    kth_distances = np.linalg.norm(kth_neighbors - positions, axis=1)
    return kth_distances


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


# Cache the kernels to avoid recomputation and significantly improve runtime
@lru_cache(maxsize=256)
def create_gaussian_kernel(sigma):
    x = np.arange(-3 * sigma, 3 * sigma + 1, 1)
    y = np.arange(-3 * sigma, 3 * sigma + 1, 1)
    xx_kernel, yy_kernel = np.meshgrid(x, y)
    kernel = np.exp(-(xx_kernel**2 + yy_kernel**2) / (2.0 * sigma**2)) / (2.0 * np.pi * sigma**2)
    return kernel


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
                A [npx, npx] image will be produced later.
                It accept 'auto' parameter, which will automatically decide the pixel number to include all particles.
    AR      : Angular resolution. Type: arcsec. Default: None
                If AR is None and npx is auto, we will force the npx to be 512
                If AR is None, it will automatically recalculated in the code for output.
                At z=0, no matter AR is set or not, AR is always assume at z = 0.05.
    redshift: The redshift of the object at. Default: None.
                If None, redshift from simulation data will be used.
                This will be moved to 0.10 if simulation redshift is used and equal to 0.
                Note this redshift does not change anything of the simulation physical particle positions, only shift the object to this redshift for observing.
    zthick  : The thickness in projection direction. Default: None.
                If None, use all data from cutting region. Otherwise set a value in simulation
                length unit (kpc/h normally), then a slice of data [center-zthick, center+zthick]
                will be used to make the y-map.
    SP      : Faked sky positions in [RA (longitude), DEC (latitude)] in degrees.
                Default: [194.95, 27.98], Coma' position.
                If [x,y,z] (len(SP) == 3) of the Earth position in the pos coordinate is given,
                The pos - [x,y,z] are taken as the J2000 3D coordinates and converted into RA, DEC.
    unit    : is the input wdata in luminosity, flux or magnitude? Default: flux, assume in ab mag.
                Set this to true if particles' luminosity are given in wdata.
    ksmooth : An integer representing the k in kNN Gaussian smoothing. 1 sigma for the Gaussian is set to the distance to the kth neighbour in 3D space.
                Recommended value: somewhere between 20 and 80.
                If k>0, you set the smoothing length of a particle to be the distance between the particle and its kth nearest neighbour.
                If k=0, pymgal does not perform smoothing.
                If k<0, throw an error.
    lsmooth:  An array of floats where each float is the smoothing length (ie kNN distance) for a given particle.
                Default: None, but can be set to a precomputed array to avoid redundant calculations
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
    Pdata = pymgal.projection(part_lum, simu_data, npx=1024, unit='flux')
    Pdata.write_fits_image("filename.fits")
    """

    def __init__(self, data, simd, axis="z", npx=512, AR=None, redshift=None, zthick=None,
                 SP=[194.95, 27.98], unit='flux', ksmooth=0, lsmooth=None, outmas=False, outage=False, outmet=False):

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
        if len(SP) == 2:
            self.sp = SP
            self.cc = simd.center
        elif len(SP) == 3:
            self.sp = False
            self.cc = SP
        else:
            raise ValueError("SP length should be either 2 or 3!")
        self.rr = simd.radius/simd.cosmology.h / (1.+ simd.redshift)   # to physical in simulation time
        self.flux = unit
        self.omas = outmas
        self.oage = outage
        self.omet = outmet
        self.ksmooth = ksmooth
        if ksmooth < 0:
            raise ValueError("ksmooth should be a value between 0 and 100 inclusively")
        
        self.lsmooth = lsmooth
        self.zthick = zthick
        if zthick is not None:
            self.zthick /= (simd.cosmology.h / (1.+ simd.redshift))
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

        pos = np.copy(s.S_pos) / s.cosmology.h / (1.+ s.redshift)  # to assumed physical
        pos -= self.cc / s.cosmology.h / (1.+ s.redshift)
        center = s.center / s.cosmology.h / (1.+ s.redshift)
        kth_distance = 0  # initialize the variable and update it if self.ksmooth > 0

        # Calculate the smoothing lengths array if it is not precomputed
        if self.ksmooth > 0 and self.lsmooth is None:
            print("Computing kNN distances")
            self.lsmooth = knn_distance(pos, self.ksmooth)


        if isinstance(self.axis, type('')):
            if self.axis.lower() == 'y':  # x-z plane
                pos = pos[:, [0, 2, 1]]
                if self.sp is False:
                    self.cc[0], self.cc[1], self.cc[2] = center[0] - self.cc[0],\
                        center[2] - self.cc[2], center[1] - self.cc[1]
            elif self.axis.lower() == 'x':  # y - z plane
                pos = pos[:, [1, 2, 0]]
                if self.sp is False:
                    self.cc[0], self.cc[1], self.cc[2] = center[1] - self.cc[1],\
                        center[2] - self.cc[2], center[0] - self.cc[0]
            else:
                if self.axis.lower() != 'z':  # project to xy plane
                    raise ValueError("Do not accept this value %s for projection" % self.axis)
                if self.sp is False:
                    self.cc = center - self.cc
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
                    self.cc = np.dot(center - self.cc, Rxyz)
            else:
                raise ValueError(
                    "Do not accept this value %s for projection" % self.axis)
        elif isinstance(self.axis, type(np.array([]))):
            if len(self.axis.shape) == 1:  # only vector from the line of sight.
                normed_axis = self.axis/np.sqrt(np.sum(self.axis**2))
                pos = np.cross(normed_axis, np.cross(pos, normed_axis))
                if self.sp is False:
                    self.cc = np.cross(normed_axis, np.cross(center - self.cc, normed_axis))
            elif len(self.axis.shape) == 2:
                if self.axis.shape[0] == self.axis.shape[1] == 3:
                    pos = np.dot(pos, self.axis)
                    if self.sp is False:
                        self.cc = np.dot(center - self.cc, self.axis)
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

        if self.sp is not False:  # noraml projection
            minx = -(self.npx + 1) * self.pxsize / 2
            maxx = +(self.npx + 1) * self.pxsize / 2
            miny = -(self.npx + 1) * self.pxsize / 2
            maxy = +(self.npx + 1) * self.pxsize / 2
            xx = np.arange(minx, maxx, self.pxsize)
            yy = np.arange(miny, maxy, self.pxsize)
        else:  # we do real projection by transferring into RA, Dec
            SC = SkyCoord(pos[:, 0], pos[:, 1], pos[:, 2], unit='kpc',
                          representation='cartesian')
            SC = SC.transform_to('icrs')
            pos[:, 0], pos[:, 1] = SC.ra.degree, SC.dec.degree

            SC = SkyCoord(self.cc[0], self.cc[1], self.cc[2], unit='kpc',
                          representation='cartesian')
            SC = SC.transform_to('icrs')
            self.sp = [SC.ra.degree, SC.dec.degree]

            minx = self.sp[0] - (self.npx + 1) * self.ar / 2.
            maxx = self.sp[0] + (self.npx + 1) * self.ar / 2.
            miny = self.sp[1] - (self.npx + 1) * self.ar / 2.
            maxy = self.sp[1] + (self.npx + 1) * self.ar / 2.
            xx = np.arange(minx, maxx, self.ar)
            yy = np.arange(miny, maxy, self.ar)


        x_bins = np.digitize(pos[:, 0], xx)                                                     # define x and y bins
        y_bins = np.digitize(pos[:, 1], yy)

        self.cc = center  # real center in the data
        L_sun = 3.826e33 # convert from solar luminosities to erg/s
        # If smoothing is set to off
        if self.ksmooth == 0:
            for i in d.keys():
                if self.flux.lower() == 'luminosity':  # luminosity
                    self.outd[i] = np.histogram2d(pos[:, 0], pos[:, 1], bins=[xx, yy], weights=d[i]*L_sun)[0]
                elif self.flux.lower() == 'flux': # flux
                    self.outd[i] = np.histogram2d(pos[:, 0], pos[:, 1], bins=[xx, yy], weights=d[i]*L_sun*100./s.cosmology.luminosity_distance(self.z).to('pc').value**2)[0]
                else: # ab mag
                    self.outd[i] = np.ma.log10(np.histogram2d(pos[:, 0], pos[:, 1], bins=[xx, yy], weights=10**(d[i]/-2.5))[0])
                    self.outd[i] = -2.5*self.outd[i].filled(self.outd[i].min()/2.)

        # If smoothing is set to on and we have at least one filter to smooth
        elif (self.ksmooth > 0) and (len(d.keys()) > 0):  
            sample_hist =  np.histogram2d(pos[:, 0], pos[:, 1], bins=[xx, yy], weights=list(d.values())[0])[0] # define a sample hist using the weights from whichever filter comes first
            smoothed_hists = None # initialize the set of smoothed histograms 
            if self.flux.lower() == 'luminosity':  # luminosity
                 d_lum = {key: d[key]*L_sun for key in d}
                 smoothed_hists = gaussian_smoothing(self.lsmooth, d_lum, sample_hist, x_bins, y_bins, self.pxsize)
            elif self.flux.lower() == 'flux': # flux
                d_flux = {key: d[key]*L_sun*100./s.cosmology.luminosity_distance(self.z).to('pc').value**2 for key in d}
                smoothed_hists = gaussian_smoothing(self.lsmooth, d_flux, sample_hist, x_bins, y_bins, self.pxsize)
            else: # ab mag
                d_mag = {key: 10**(d[key]/-2.5) for key in d}
                smoothed_hists = gaussian_smoothing(self.lsmooth, d_mag, sample_hist, x_bins, y_bins, self.pxsize)

            # Set the histograms to match the filters
            for i in d.keys():
                if self.flux.lower() == 'magnitude':
                    self.outd[i] = np.ma.log10(smoothed_hists[i])
                    self.outd[i] = -2.5*self.outd[i].filled(self.outd[i].min()/2.)
                else:
                    self.outd[i] = smoothed_hists[i]


        

        # Now grid the data
        # pmax, pmin = np.max(self.S_pos, axis=0), np.min(self.S_pos, axis=0)
        # grid_x, grid_y = np.mgrid[pmin[0]:pmax[0]:nx, pmin[1]:pmax[1]:nx]
        if self.omas or self.oage or self.omet:
            #unsmooth_hist = np.histogram2d(pos[:, 0], pos[:, 1], bins=[xx, yy], weights=s.S_mass)[0]
            
            if self.ksmooth > 0:
                max_sigma = 5/self.pxsize # Set the max standard deviation for the smoothing gaussian to be no more than the gravitational softening length of 5 kpc/h
                sample_hist = np.histogram2d(pos[:, 0], pos[:, 1], bins=[xx, yy], weights=s.S_mass)[0]
                mass_hist_dict = gaussian_smoothing(self.lsmooth, {"Mass": s.S_mass}, sample_hist, x_bins, y_bins, self.pxsize, max_kernel=max_sigma)
                self.outd["Mass"] = mass_hist_dict["Mass"]

            else:
                self.outd["Mass"] = np.histogram2d(pos[:, 0], pos[:, 1], bins=[xx, yy], weights=s.S_mass)[0]
            
        if self.oage:
            unsmooth_mass_hist = np.histogram2d(pos[:, 0], pos[:, 1], bins=[xx, yy], weights=s.S_mass)[0]
            ids = unsmooth_mass_hist > 0
            self.outd["Age"] = np.histogram2d(pos[:, 0], pos[:, 1], bins=[xx, yy], weights=s.S_mass * s.S_age)[0]
            self.outd["Age"][ids] /= unsmooth_mass_hist[ids]
                
        if self.omet:
            unsmooth_mass_hist = np.histogram2d(pos[:, 0], pos[:, 1], bins=[xx, yy], weights=s.S_mass)[0]
            ids = unsmooth_mass_hist > 0
            self.outd["Metal"] = np.histogram2d(pos[:, 0], pos[:, 1], bins=[xx, yy], weights=s.S_mass * s.S_metal)[0]
            self.outd["Metal"][ids] /= unsmooth_mass_hist[ids]

        

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
            hdu.header["PIXVAL"] = self.flux
            hdu.header.comments["PIXVAL"] = 'in flux[ergs/s/cm^2], lumi[ergs/s] or mag.'
            hdu.header["ORAD"] = float(self.rr)
            hdu.header.comments["ORAD"] = 'Rcut in physical for the image.'
            hdu.header["REDSHIFT"] = float(self.z)
            hdu.header.comments["REDSHIFT"] = 'The redshift of the object being put to'
            hdu.header["PSIZE"] = float(self.pxsize)
            hdu.header.comments["PSIZE"] = 'The pixel size in physical at simulation time'
            hdu.header["AGLRES"] = float(self.ar*3600.)
            hdu.header.comments["AGLRES"] = '\'observation\' angular resolution in arcsec'
            hdu.header["ORIGIN"] = 'PymGal'
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


