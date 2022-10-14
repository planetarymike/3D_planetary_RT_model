import numpy as np
import datetime
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.io import fits
import idl_colorbars
import os

import matplotlib as mpl
import matplotlib.pyplot as plt

import matplotlib.font_manager as fm
mpl.rcParams['font.sans-serif'] = 'Gillius ADF'
mpl.rcParams['font.family'] = "sans-serif"

emus_files_path = '/home/mike/Documents/EMM/data/emus' # where the EMUS files live
quicklook_save_path = '/home/mike/Documents/EMM/quicklooks' # where to save generated quicklooks
quicklook_code_location = '/home/mike/Documents/EMM/emus_sci/utils/quicklook/' # where the quicklook.py file is

current_directory = os.getcwd()
os.chdir(quicklook_code_location)
import quicklook
os.chdir(current_directory)

import datetime
from scipy.io.idl import readsav
import spiceypy as spice

# routine to get EUVM Lyman alpha line center brightness
euvm_l2b_filename = '/home/mike/Documents/MAVEN/EUVM/mvn_euv_l2b_orbit_merged_v14_r02.sav'

def get_solar_lyman_alpha(input_datetime):
    """Compute the EUVM-measured solar Lyman alpha line center values
    for the input et_list. Uses orbit-averaged EUVM l2b file synced
    from MAVEN SDC. Requires spicepy (but no SPICE kernels) to convert
    ET to datetime.

    Parameters
    ----------
    input_et : float
        SPICE ETs at which solar lyman alpha brightness is needed.

    Returns
    -------
    lya_interp : float or numpy.nan
        Interpolated MAVEN EUVM line center Lyman alpha brightness value
        at Mars. Units are photons/cm2/s/nm. If no EUVM data are 
        available within 2 days of the requested ET, returns np.nan.

    """

    # load the EUVM data
    euvm_l2b = readsav(euvm_l2b_filename)['mvn_euv_l2_orbit'].item()
    # start time of EUVM orbit
    euvm_datetime = [datetime.datetime.fromtimestamp(t) for t in euvm_l2b[0]]

    euvm_lya = euvm_l2b[2][2] # median orbit EUVM brightness
    euvm_mars_sun_dist = euvm_l2b[5]

    # we need to remove the timezone info to compare with EUVM times
    input_datetime = input_datetime.replace(tzinfo=None)

    # interpolate the EUVM data if it's close enough in time
    euvm_idx = np.searchsorted(euvm_datetime, input_datetime) - 1
    
    if euvm_idx == len(euvm_datetime)-1:
        raise ValueError("End time of EUVM data "
                         f"{datetime.datetime.isoformat(euvm_datetime[-1])} "
                         "is before requested time "
                         f"{datetime.datetime.isoformat(input_datetime)}.")
    
    days_adjacent = 2
    erly_cutoff = input_datetime - datetime.timedelta(days=days_adjacent)
    late_cutoff = input_datetime + datetime.timedelta(days=days_adjacent)

    after_erly_cutoff = euvm_datetime[euvm_idx  ] > erly_cutoff
    befor_late_cutoff = euvm_datetime[euvm_idx+1] < late_cutoff
    if after_erly_cutoff and befor_late_cutoff:
        inpt_timediff = input_datetime - euvm_datetime[euvm_idx]
        euvm_timediff = euvm_datetime[euvm_idx+1] - euvm_datetime[euvm_idx]

        interp_frac = inpt_timediff / euvm_timediff

        lya_interp  = (interp_frac*(euvm_lya[euvm_idx + 1]
                                    - euvm_lya[euvm_idx])
                       + euvm_lya[euvm_idx])
        dist_interp = (interp_frac*(euvm_mars_sun_dist[euvm_idx + 1]
                                    - euvm_mars_sun_dist[euvm_idx])
                       + euvm_mars_sun_dist[euvm_idx])
        dist_interp = dist_interp / 1.496e8  # convert dist_interp to AU

        # new we have the band-integrated value measured at Mars, we
        # need to convert back to Earth, then get the line center flux
        # using the power law relation of Emerich+2005
        lya_interp *= dist_interp**2
        # this is now in W/m2 at Earth. We need to convert to ph/cm2/s
        phenergy = 1.98e-25/(121.6e-9)  # energy of a lyman alpha photon in J
        lya_interp /= phenergy
        lya_interp /= 1e4  # convert to /cm2
        # we're now in ph/cm2/s

        # Use the power law relation of Emerich+2005:
        lya_interp = 0.64*((lya_interp/1e11)**1.21)
        lya_interp *= 1e12
        # we're now in ph/cm2/s/nm

        # convert back to Mars
        lya_interp /= dist_interp**2
    else:
        lya_interp = np.nan

    return lya_interp

def solar_flux_to_g_factor(solar_flux_line_center, wavelength=121.56, f_value=0.416):
    # converts EUVM fluxes in ph/cm2/s/nm to a g factor
    wavelength = wavelength*1e-7 # cm
    cross_section_total = 2.647e-2 * f_value  # cm^2 Hz
    clight = 3e10  # cm/s
    g_factor = (solar_flux_line_center * 1e7 # nm / cm
                * wavelength * wavelength) / clight
    # now we're in photons/s/cm2/Hz
    g_factor *= cross_section_total
    return g_factor

all_os2 = sorted(quicklook.find_files(data_directory=emus_files_path, mode='os2', level='l2a'))
all_os2_datetime = np.array([datetime.datetime.strptime(f.split("_")[3], "%Y%m%dt%H%M%S") for f in all_os2])
all_os2_timegap = np.diff(all_os2_datetime)
all_os2_timegap = np.array([t.total_seconds() for t in all_os2_timegap])
os2_groups = np.split(all_os2, np.where(all_os2_timegap>1e5)[0]+1)[1:]

from scipy.spatial.transform import Rotation as R
def convert_sun_mme_to_mars_ecliptic_AU(sun_mme_vecs):
    j2000_to_mme_mat = np.array([[ 0.6732521982472339,       0.7394129276360180,       0.0000000000000000],
                                 [-0.5896387605430040,       0.5368794307891331,       0.6033958972853946],
                                 [ 0.4461587269353556,      -0.4062376142607541,       0.7974417791532832]])
    j2000_to_eclip_mat = R.from_rotvec(np.pi/180. * 23.44 * np.array([1, 0, 0])).as_matrix()
    
    # convert to J2000
    outvecs = np.dot(np.linalg.inv(j2000_to_mme_mat), -np.transpose(sun_mme_vecs))
    outvecs = np.dot(j2000_to_eclip_mat, outvecs)
    
    return np.transpose(outvecs)/1.496e8

class EMUS_H_corona_data:
# a class to extract data from a given FITS file
    
    # filename
    fname = ""
    
    # dimensions of this file (n_int, n_spa)
    shape = np.array([])
    
    # integration ET (n_int)
    et = np.array([])
    utc = np.array([])
    datetime = np.array([])
        
    # mean time and solar brightness values at this time
    mean_datetime = None
    solar_lya_flux = None
    solar_lya_g = None
    solar_lyb_g = None
    
    # MSO spacecraft position (n_int, 3)
    sc_pos = np.array([]) # km
    
    # lines of sight (n_int, n_spa)
    los = np.array([])
    
    # minimum ray height altitude
    mrh = np.array([])
    
    # variables for IPH calculation
    # ra and dec (n_int, n_spa)
    ra = np.array([])
    dec = np.array([])
    # mars ecliptic coordinates ()
    mars_eclip = np.array([])
    
    # observed brightness (n_int, n_spa)
    lya_R = np.array([])
    lya_counts = np.array([])
    lyb_R = np.array([])
    lyb_counts = np.array([])
    
    # observed brightness uncertainty
    lya_R_unc = np.array([])
    lya_counts_unc = np.array([])
    lyb_R_unc = np.array([])
    lyb_counts_unc = np.array([])
    
    # conversion between R and counts 
    lya_counts_per_R = np.array([])
    lyb_counts_per_R = np.array([])
    
    # mask for good obs to fit
    mask = np.array([])
    
    def __init__(self, filename):
        self.fname = filename
        myfits = fits.open(filename)
        
        self.shape = np.array(myfits['CAL'].data['RADIANCE'].shape[:2])
        
        # time of each integration
        self.et = myfits['TIME'].data['TIME_ET']
        self.utc = myfits['TIME'].data['TIME_UTC']
        self.datetime = np.array([datetime.datetime.fromisoformat(t) 
                                  for t in self.utc])
        self.datetime_mean = \
          datetime.datetime.fromtimestamp(np.mean([t.timestamp() 
                                                   for t in self.datetime]))
        
        # EUVM solar Lyman alpha line center brightness for this file
        self.solar_lya_flux = get_solar_lyman_alpha(self.datetime_mean)
        self.solar_lya_g = solar_flux_to_g_factor(self.solar_lya_flux,
                                                   wavelength=121.56,
                                                   f_value=0.41641)
        self.solar_lyb_g = solar_flux_to_g_factor(self.solar_lya_flux/66.,
                                                   wavelength=102.57,
                                                   f_value=0.079142)
        
        # spacecraft geometry
        self.sc_pos = myfits['SC_GEOM'].data['V_SC_POS_MSO']
        
        # LOS MSO vector
        self.los = np.transpose(myfits['FOV_GEOM'].data['VEC_MSO'][:,:,4,:],
                                (0,2,1))
        
        # MRH altitude
        self.mrh = myfits['FOV_GEOM'].data['MRH_ALT'][:,4,:]
        
        # coordinates for IPH calculation
        self.ra  = myfits['FOV_GEOM'].data['RA'][:,4,:]
        self.dec = myfits['FOV_GEOM'].data['DEC'][:,4,:]
        self.mars_eclip = myfits['SC_GEOM'].data['V_SUN_POS_MME']
        self.mars_eclip = \
          convert_sun_mme_to_mars_ecliptic_AU(self.mars_eclip)
        
        # extract Lyman alpha and Lyman beta brightness in R and counts
        self.lya_R, self.lya_R_unc = \
          quicklook.integrate_swath_radiance(myfits, 
                                             [121.6 - 2.0, 121.6 + 2.0],
                                             correct_flatfield=True,
                                             return_unc=True)
        
        self.lya_counts, self.lya_counts_unc = \
          quicklook.integrate_swath_radiance(myfits, 
                                             [121.6 - 2.0, 121.6 + 2.0],
                                             correct_flatfield=True,
                                             return_radiance=False,
                                             return_unc=True)
        
        self.lyb_R, self.lyb_R_unc = \
          quicklook.integrate_swath_radiance(myfits, 
                                             [102.6 - 2.0, 102.6 + 2.0],
                                             correct_flatfield=True,
                                             return_unc=True)
        
        self.lyb_counts, self.lyb_counts_unc = \
          quicklook.integrate_swath_radiance(myfits, 
                                             [102.6 - 2.0, 102.6 + 2.0],
                                             correct_flatfield=True,
                                             return_radiance=False,
                                             return_unc=True)
        
        # counts >135.6 nm and between Lyman alpha and beta,
        # to identify off-disk stars
        self.star_counts = \
          quicklook.integrate_swath_radiance(myfits, 
                                             [135.6 - 2.0, 160.],
                                             correct_flatfield=True,
                                             return_radiance=False)
        self.star_counts += \
          quicklook.integrate_swath_radiance(myfits, 
                                             [102.5 + 2.0, 121.6 - 2.0],
                                             correct_flatfield=True,
                                             return_radiance=False)
        
        self.remove_off_sky()
        self.compute_mask()
        
        # determine conversion between R and counts for this file
        with np.errstate(divide='ignore', invalid='ignore'):
            self.lya_counts_per_R = np.nanmedian(self.lya_counts/self.lya_R)
            self.lyb_counts_per_R = np.nanmedian(self.lyb_counts/self.lyb_R)

    def remove_off_sky(self):
        sky_slit_pos = np.where(np.logical_not(np.all(np.isnan(self.ra),
                                                      axis=0)))[0]
        self.fits_sky_slit_pos = sky_slit_pos
        self.shape[1] = len(sky_slit_pos)
        
        self.los = self.los[:, sky_slit_pos]
        self.mrh = self.mrh[:, sky_slit_pos]
        self.ra = self.ra[:, sky_slit_pos]
        self.dec = self.dec[:, sky_slit_pos]
        
        self.lya_R = self.lya_R[:, sky_slit_pos]
        self.lya_counts = self.lya_counts[:, sky_slit_pos]
        self.lya_R_unc = self.lya_R_unc[:, sky_slit_pos]
        self.lya_counts_unc = self.lya_counts_unc[:, sky_slit_pos]
        
        self.lyb_R = self.lyb_R[:, sky_slit_pos]
        self.lyb_counts = self.lyb_counts[:, sky_slit_pos]       
        self.lyb_R_unc = self.lyb_R_unc[:, sky_slit_pos]
        self.lyb_counts_unc = self.lyb_counts_unc[:, sky_slit_pos]
        
        self.star_counts = self.star_counts[:, sky_slit_pos]
        
    def compute_mask(self):
        minimum_altitude = 500 # km
        star_counts_max = 50 # counts
        
        self.mask = np.logical_and(self.mrh > minimum_altitude,
                                   self.star_counts < star_counts_max)

import scipy.stats
import time

os.chdir('/home/mike/Documents/Mars/3D_planetary_RT_model/python/')
from py_corona_sim_cpu import Pyobservation_fit
os.chdir(current_directory)

# a class that stores all model input (EMM+Lyman alpha), and performs a fit
class EMUS_H_corona_fit:   
    def __init__(self, filename_list):
        self.file_data = []
        self.file_shape = []
        self.file_mask = []
        self.file_n_los = np.array([])
        self.file_startend = np.array((0,2))

        # compiled list of data to fit from all OSs
        # reshaped so there's one entry per LOS in all files
        self.obs_solar_lya_g = np.array([])
        self.obs_solar_lyb_g = np.array([])

        self.obs_sc_pos = np.empty((0,3)) # km
        self.obs_los = np.empty((0,3))

        self.obs_ra = np.array([])
        self.obs_dec = np.array([])
        self.obs_mars_eclip = np.empty((0,3))

        # observed brightness (n_int, n_spa)
        self.obs_lya_counts = np.array([])
        self.obs_lyb_counts = np.array([])

        # conversion between R and counts 
        self.lya_counts_per_R = np.array([])
        self.lyb_counts_per_R = np.array([])

        # model object
        self.corona_model = Pyobservation_fit()
        self.iph_model_init = False
        
        for f in filename_list:
            self.__load_file_data(f)
        
        # compute start and end indices for the files
        file_los_end = np.int64(np.cumsum(self.file_n_los))
        file_los_start = np.concatenate([[0], file_los_end[:-1]])
        self.file_startend = np.transpose([file_los_start,
                                           file_los_end])
        
        self._model_obs_setup()
        
    def __load_file_data(self, filename):
        # load data for a single input file
        tfile = EMUS_H_corona_data(filename)
        self.file_data.append(tfile)

        tfile_shape = tfile.shape
        self.file_shape.append(tfile_shape)

        tfile_mask = tfile.mask
        self.file_mask.append(tfile_mask)

        n_tfile_los = int(np.sum(tfile_mask))
        self.file_n_los = np.append(self.file_n_los,
                                    n_tfile_los)

        tfile_solar_lya_g = np.repeat(tfile.solar_lya_g, 
                                          n_tfile_los)
        self.obs_solar_lya_g = np.append(self.obs_solar_lya_g,
                                         tfile_solar_lya_g)
        tfile_solar_lyb_g = np.repeat(tfile.solar_lyb_g, n_tfile_los)
        self.obs_solar_lyb_g = np.append(self.obs_solar_lyb_g,
                                         tfile_solar_lyb_g)

        tfile_sc_pos = np.repeat(tfile.sc_pos[:, np.newaxis, :],
                                tfile_shape[1],
                                1)
        self.obs_sc_pos = np.append(self.obs_sc_pos,
                                    tfile_sc_pos[tfile_mask], 
                                    axis=0)
        self.obs_los = np.append(self.obs_los,
                                 tfile.los[tfile_mask], 
                                 axis=0)

        self.obs_ra = np.append(self.obs_ra,
                                 tfile.ra[tfile_mask])
        self.obs_dec = np.append(self.obs_dec,
                                 tfile.dec[tfile_mask])
        tfile_mars_eclip = np.repeat(tfile.mars_eclip[:, np.newaxis, :],
                                    tfile_shape[1],
                                    1)
        self.obs_mars_eclip = np.append(self.obs_mars_eclip,
                                        tfile_mars_eclip[tfile_mask], 
                                        axis=0)

        self.obs_lya_counts = np.append(self.obs_lya_counts,
                                        tfile.lya_counts[tfile_mask])
        self.obs_lyb_counts = np.append(self.obs_lyb_counts,
                                        tfile.lyb_counts[tfile_mask])

        tfile_lya_counts_per_R = np.repeat(tfile.lya_counts_per_R, n_tfile_los)
        self.lya_counts_per_R = np.append(self.lya_counts_per_R,
                                          tfile_lya_counts_per_R)
        tfile_lyb_counts_per_R = np.repeat(tfile.lyb_counts_per_R, n_tfile_los)
        self.lyb_counts_per_R = np.append(self.lyb_counts_per_R,
                                          tfile_lyb_counts_per_R)
     
    def _file_split(self, q):
        # split input data member into lists for each file
        return [q[i_low:i_high] for i_low, i_high in self.file_startend]
    
    def _file_split_and_shape(self, q, pad=np.nan):
        # split input data member into arrays of the same shape as the original files, padding with 'pad' for excluded values
        split = self._file_split(q)
        retarr = [np.empty(s) for s in self.file_shape]
        for i in range(len(retarr)):
            retarr[i][self.file_mask[i]] = split[i]
            retarr[i][np.logical_not(self.file_mask[i])] = pad
        
        return retarr
    
    def _model_obs_setup(self):
        # load observation geometry into model
        self.corona_model.add_observation(self.obs_sc_pos*1e5, self.obs_los)
        time.sleep(0.1)        
        self.corona_model.set_g_factor([np.mean(self.obs_solar_lya_g),
                                        np.mean(self.obs_solar_lyb_g)])
        time.sleep(0.1)
        self._model_obs_iph()
        
    def _model_obs_iph(self):
        # routine to obtain IPH brightness for RA/Dec pointing
        if not self.iph_model_init:
            # the Quemerais IPH model is slow, call it only once
            print('modeling IPH, this might take a few minutes...')
            time.sleep(0.1)
            self.corona_model.add_observation_ra_dec(np.mean(self.obs_mars_eclip, axis=0), 
                                                     self.obs_ra, 
                                                     self.obs_dec)
            time.sleep(0.1)
            print('...done modeling IPH.')
        
    def model_obs_R(self,
                    nHavg, Tavg,
                    model_scale_factor_lya,
                    model_scale_factor_lyb):
        # spherically symmetric model
        self.corona_model.generate_source_function(nHavg, Tavg)
        
        # get brightness output
        model_R = self.corona_model.brightness()
        model_lya_R = (1e3 # R / kR
                       *model_scale_factor_lya
                       *self.obs_solar_lya_g/np.mean(self.obs_solar_lya_g)
                       *model_R[0])
        model_lyb_R = (1e3 # R / kR
                       *model_scale_factor_lyb
                       *self.obs_solar_lyb_g/np.mean(self.obs_solar_lyb_g)
                       *model_R[1])
        
        return [model_lya_R, model_lyb_R]
        
    def model_obs_counts(self, 
                         nHavg, Tavg,
                         model_scale_factor_lya,
                         model_scale_factor_lyb):
        
        model_lya_R, model_lyb_R = self.model_obs_R(nHavg, Tavg,
                                                    model_scale_factor_lya,
                                                    model_scale_factor_lyb)
        
        model_lya_counts = model_lya_R*self.lya_counts_per_R
        model_lyb_counts = model_lyb_R*self.lyb_counts_per_R
        
        return (model_lya_counts, model_lyb_counts)
     
    def neg_log_likelihood(self, model_params):
        nHavg = model_params[0]
        Tavg = model_params[1]
        model_scale_factor_lya = model_params[-2]
        model_scale_factor_lyb = model_params[-1]
        
        model_lya_counts, model_lyb_counts = \
            self.model_obs_counts(nHavg, Tavg,
                                  model_scale_factor_lya,
                                  model_scale_factor_lyb)
        
        # poisson likelihood for counts
        def poisson_logl(observed, expected):
            logl = (observed*np.log(expected)
                    - expected
                    - scipy.special.gammaln(observed+1))
            return logl
        
        neg_logl = -np.nansum(poisson_logl(self.obs_lya_counts, model_lya_counts))
        neg_logl += -np.nansum(poisson_logl(self.obs_lyb_counts, model_lyb_counts))
        
        if np.isnan(neg_logl):
            warnings.warn(f"logl is NaN at inputs {model_params}")
            neg_logl = np.inf
        
        # Gaussian prior for model scale factors
        neg_logl += (model_scale_factor_lya/0.2)**2
        neg_logl += (model_scale_factor_lyb/0.2)**2
        
        return neg_logl
        
    def model_fit_nelder_mead(self):
        from scipy.optimize import minimize, Bounds

        param_bounds = [[1.0e4, 1.0e7],
                        [100., 400.],
                        [0.5, 1.5],
                        [0.25, 4.0]]

        # do the minimization
        bestfit = minimize(self.neg_log_likelihood,
                           (1.0e5, 300., 1.0, 1.0),
                           method='Nelder-Mead',
                           bounds=param_bounds)
        
        return bestfit

test_emus_fit = EMUS_H_corona_fit(os2_groups[2][:1])
