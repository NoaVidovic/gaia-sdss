import gaiaxpy as gxp
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import sdss

from astroML.sum_of_norms import sum_of_norms, norm
from astroquery.sdss import SDSS as aq_sdss
from scipy.optimize import curve_fit
from sys import argv



# UTIL FUNCTIONS
def get_script_path():
    return os.path.dirname(os.path.realpath(argv[0]))

def get_q():
    data_q = os.path.join(get_script_path(), '../data/q_median.npy')
    src_q  = os.path.join(get_script_path(), './q_median.npy')

    if os.path.exists(data_q):
        return np.load(data_q)
    if os.path.exists(src_q):
        return np.load(src_q)
    return None


def get_main_table():
    data_mt = os.path.join(get_script_path(), '../data/main_table.csv')
    src_mt  = os.path.join(get_script_path(), './main_table.csv')

    if os.path.exists(data_mt):
        return pd.read_csv(data_mt)
    if os.path.exists(src_mt):
        return pd.read_csv(src_mt)
    return None


def get_figure_path():
    data_fs = os.path.join(get_script_path(), '../figures')
    src_fs  = os.path.join(get_script_path(), './figures')

    if os.path.exists(data_fs):
        return data_fs

    if not os.path.exists(src_fs):
        os.makedirs(src_fs)

    return src_fs



## Delta lambda as a function of lambda for Gaia spectrum
def deltaXP(wvl):
    y1 = lambda x: 3.30 * (x/330)**2.45
    y2 = lambda x: 6.40 * (x/640)**1.72
    y_mid = lambda x: (y2(680) * (x - 640) + y1(640) * (680 - x))/40
    
    return np.piecewise(wvl, \
                        [wvl < 640, (640 <= wvl) & (wvl <= 680), wvl > 680], \
                        [y1, y_mid, y2])


# CONSTANTS

DF = get_main_table()
Q_MEDIAN = get_q()
FIGURE_PATH = get_figure_path()
MIN_WVL, MAX_WVL = 390, 900
NUM_NORMS = 800

## plot settings
SMALL_FONT_SIZE = 16
MEDIUM_FONT_SIZE = 20
LARGE_FONT_SIZE = 20

plt.rc('font', size=SMALL_FONT_SIZE)          # default text sizes
plt.rc('axes', titlesize=LARGE_FONT_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_FONT_SIZE)    # fontsize of the x and y labels
plt.rc('axes', linewidth=2)                   # linewidth of the graph borders
plt.rc('xtick', labelsize=SMALL_FONT_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_FONT_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_FONT_SIZE)    # legend fontsize
plt.rc('figure', titlesize=LARGE_FONT_SIZE)   # fontsize of the figure title

SDSS_RAW_COLOR = '#AAAAAA'
SDSS_FIT_COLOR = '#FF4444'
GAIA_COLOR = '#4444FF'

# main function

def plot_spectra(*, gaia_id=None, sdss_id=None, savefig=False):
    if gaia_id is None and sdss_id is None:
        print('No ID provided')
        return
    elif sdss_id is None:    
        try:
            sdss_id = DF.loc[DF['source_id'] == gaia_id]['specObjId'].to_list()[0]
        except:
            print('Could not find an SDSS ID corresponding to the provided Gaia ID')
            return
    elif gaia_id is None:
        try:
            gaia_id = DF.loc[DF['specObjId'] == sdss_id]['source_id'].to_list()[0]
        except:
            print('Could not find an Gaia ID corresponding to the provided SDSS ID')
            return
            
    # get sdss data
    try:
        sp = sdss.SpecObj(int(sdss_id))
        data = aq_sdss.get_spectra(plate=sp.plate, mjd=sp.mjd, fiberID=sp.fiberID)
    except:
        print('Could not get SDSS data')
        return
    
    spec_data = data[0][1].data

    sdss_sampling = 10 ** spec_data['loglam'] / 10  # Convert log wavelength to linear and Å to nm
    sdss_flux = spec_data['flux'] * 1e-19 # Convert SDSS units (1e-17 in cgs) to Gaia units (SI)
    
    # get gaia data and calibrate using sdss sampling
    gaia_flux = gxp.calibrate([gaia_id], sampling=sdss_sampling, truncation=True)[0]['flux'][0]

    weights, fit_rms, locs, widths = sum_of_norms(sdss_sampling, sdss_flux, NUM_NORMS,
                                                 spacing='linear',
                                                 full_output=True)
    
    # convolve sdss with sigma_gaia
    k = 0.474  # obtained experimentally
    sigma_conv = np.sqrt(widths**2 + k**2 * deltaXP(sdss_sampling[:, None])**2)  # convolve gaussian sigmas with sigma_gaia
    sdss_flux_fit = (weights * norm(sdss_sampling[:, None], locs, widths)).sum(1)
    sdss_conv = (weights * norm(sdss_sampling[:, None], locs, sigma_conv)).sum(1)
    
    # mask everything to the area of interest
    mask = (sdss_sampling > 390) & (sdss_sampling < 900)

    sampling = sdss_sampling[mask]
    sdss_flux = sdss_flux[mask]
    sdss_flux_fit = sdss_flux_fit[mask]
    gaia_flux = gaia_flux[mask]
    sdss_conv = sdss_conv[mask]

    floors = np.floor(sampling).astype(np.int64)
    closer_to_higher = (sampling - floors) > 0.5
    indices = floors - MIN_WVL + closer_to_higher

    gaia_corr = gaia_flux * Q_MEDIAN[indices]

    w = sdss_conv / gaia_corr
    w_norm = w / np.median(w)
    w_25 = np.quantile(w, 0.25)
    w_75 = np.quantile(w, 0.75)
    w_rrms = 0.7413 * (w_75 - w_25)

    xticks = np.arange(MAX_WVL, MIN_WVL, -20)[::-1]
    
    # plot the results
    fig, ax = plt.subplots(nrows=3, figsize=(30, 30))

    # raw data
    ax[0].plot(sampling, sdss_flux,     color=SDSS_RAW_COLOR, ls='--', lw=1, label='SDSS raw data')
    ax[0].plot(sampling, sdss_flux_fit, color=SDSS_FIT_COLOR, ls='-',  lw=2, label='SDSS Gaussian fit')
    ax[0].plot(sampling, gaia_flux,     color=GAIA_COLOR,     ls='-',  lw=2, label='Gaia raw data')

    ax[0].legend(loc=0)
    ax[0].set_xlim(MIN_WVL, MAX_WVL)
    ax[0].set_xticks(xticks)
    ax[0].set_ylabel('flux [W m$^{-2}$ nm$^{-1}$]')

    # corrected spectra
    ax[1].plot(sampling, sdss_flux, color=SDSS_RAW_COLOR, ls='--', lw=1, label='SDSS raw data')
    ax[1].plot(sampling, sdss_conv, color=SDSS_FIT_COLOR, ls='-',  lw=2, label='SDSS convolved')
    ax[1].plot(sampling, gaia_corr, color=GAIA_COLOR,     ls='-',  lw=2, label='Gaia corrected')

    ax[1].legend(loc=0)
    ax[1].set_xlim(MIN_WVL, MAX_WVL)
    ax[1].set_xticks(xticks)
    ax[1].set_ylabel('flux [W m$^{-2}$ nm$^{-1}$]')

    # w
    ax[2].axhline(1, color='#888888', ls='--', lw=2)
    ax[2].plot(sampling, w_norm, color='#000000', ls='-',  lw=2)

    ax[2].set_xlim(MIN_WVL, MAX_WVL)
    ax[2].set_xticks(xticks)
    ax[2].set_ylabel('$w$ / median($w$)')

    ax[2].text(0.97, 0.95, f"median(w) = {np.median(w):.3}", ha='right', va='top', transform=plt.gca().transAxes)
    ax[2].text(0.97, 0.85, f"rrms(w) = {w_rrms:.3}", ha='right', va='top', transform=plt.gca().transAxes)

    ax[2].set_xlabel('$\\lambda$ [nm]')
    ax[0].set_title(f"Gaia source_id {gaia_id}\n SDSS specObjId {sdss_id}")
    
    if savefig:
        plt.savefig(os.path.join(FIGURE_PATH, f'spectra_w_G{gaia_id}_S{sdss_id}.png'))

    plt.show()


if __name__ == '__main__':
    try:
        gid = int(argv[1])
    except:
        print("No ID provided")

    plot_spectra(gaia_id=gid, savefig=True)
