import gaiaxpy as gxp
import numpy as np
import os
import pandas as pd
import sdss
import warnings

from astroML.sum_of_norms import sum_of_norms, norm
from astroquery.sdss import SDSS as aq_sdss
from tqdm import tqdm

from util import get_script_path


# constants
DATA_PATH = os.path.join(get_script_path(), '..', 'data')
SPECTRA_PATH = os.path.join(DATA_PATH, 'spectra')

DF = pd.read_csv(os.path.join(DATA_PATH, 'main_table.csv'))
df_copy = DF.copy()
df_copy['sdss_flux'] = 0
df_copy['gaia_flux'] = 0
df_copy['mean_q'] = 0
df_copy['median_q'] = 0

GAIA_ID_COLNAME = 'source_id'
SDSS_ID_COLNAME = 'specObjId'
NUM_NORMS = 800


def deltaXP(wvl):
    """
    Delta lambda (lambda) for Gaia XP spectra.
    The overlap between the BP and RP portions of the spectra is joined linearly.
    """
    y1 = lambda x: 3.30 * (x/330)**2.45
    y2 = lambda x: 6.40 * (x/640)**1.72
    y_mid = lambda x: (y2(680) * (x-640) + y1(640) * (680-x)) / (680-640)

    return np.piecewise(wvl, \
                        [wvl < 640, (640 <= wvl) & (wvl <= 680), wvl > 680], \
                        [y1, y_mid, y2])


def get_q(*, gaia_id=None, sdss_id=None, k=0.474):
    """
    Gets q for a given SDSS ID.
    """
    if gaia_id is None and sdss_id is None:
        print('No ID supplied')
    elif sdss_id is None:    
        try:
            sdss_id = DF.loc[DF[GAIA_ID_COLNAME] == gaia_id][SDSS_ID_COLNAME].to_list()[0]
        except:
            print('Could not find an SDSS ID corresponding to the provided Gaia ID')
    elif gaia_id is None:
        try:
            gaia_id = DF.loc[DF[SDSS_ID_COLNAME] == sdss_id][GAIA_ID_COLNAME].to_list()[0]
        except:
            print('Could not find an Gaia ID corresponding to the provided SDSS ID')

    # SDSS spectrum
    try:
        sp = sdss.SpecObj(int(sdss_id))
        data = aq_sdss.get_spectra(plate=sp.plate, mjd=sp.mjd, fiberID=sp.fiberID)[0][1].data
    except:
        return (-1, -1, -1, -1)

    sdss_sampling = 10 ** data['loglam'] / 10
    sdss_flux = data['flux'] * 1e-19

    weights, fit_rms, locs, widths = sum_of_norms(sdss_sampling, sdss_flux, NUM_NORMS,
                                                  spacing='linear',
                                                  full_output=True)

    sigma_conv = np.sqrt(widths**2 + k**2 * deltaXP(sdss_sampling[:, None])**2)
    sdss_flux_fit = (weights * norm(sdss_sampling[:, None], locs, widths)).sum(1)
    sdss_conv = (weights * norm(sdss_sampling[:, None], locs, sigma_conv)).sum(1)

    # Gaia XP spectrum
    gaia_flux = gxp.calibrate([gaia_id], sampling=sdss_sampling, truncation=True)[0]['flux'][0]

    # integrate flux
    sdss_flux_integrated = np.trapz(sdss_flux, sdss_sampling)
    gaia_flux_integrated = np.trapz(gaia_flux, sdss_sampling)

    q = sdss_conv/gaia_flux

    if not os.path.exists(SPECTRA_PATH):
        os.makedirs(SPECTRA_PATH)

    q_path = os.path.join(SPECTRA_PATH, f'spectra_G{gaia_id}_S{sdss_id}.npy')
    np.save(q_path, np.array([sdss_sampling, sdss_flux, sdss_flux_fit, sdss_conv, gaia_flux, q]))

    return (sdss_flux_integrated, gaia_flux_integrated, np.mean(q), np.median(q))


if __name__ == '__main__':
    warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)
    for gaia_id in tqdm(DF[:][GAIA_ID_COLNAME]):
        sdss_flux, gaia_flux, mean_q, median_q = get_q(gaia_id=gaia_id)

        df_copy.loc[df_copy[GAIA_ID_COLNAME] == gaia_id, 'sdss_flux'] = sdss_flux
        df_copy.loc[df_copy[GAIA_ID_COLNAME] == gaia_id, 'gaia_flux'] = gaia_flux
        df_copy.loc[df_copy[GAIA_ID_COLNAME] == gaia_id, 'mean_q'] = mean_q
        df_copy.loc[df_copy[GAIA_ID_COLNAME] == gaia_id, 'median_q'] = median_q

    print('Data processed. Saving to main_table.csv...')
    df_copy.to_csv(os.path.join(DATA_PATH, 'main_table.csv'), index=False)
    print('Done.')

