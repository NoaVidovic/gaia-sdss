import gaiaxpy as gxp
import numpy as np
import os
import pandas as pd
import sdss
import warnings

from astroML.sum_of_norms import sum_of_norms, norm
from astroquery.sdss import SDSS as aq_sdss
from tqdm import tqdm

from util import get_script_path, decode_npy


# constants
DATA_PATH = os.path.join(get_script_path(), '..', 'data')
SPECTRA_PATH = os.path.join(DATA_PATH, 'spectra')
W_PATH = os.path.join(DATA_PATH, 'w')

DF = pd.read_csv(os.path.join(DATA_PATH, 'main_table.csv'))
df_copy = DF.copy()
df_copy['w_mean'] = 0
df_copy['w_std'] = 0
df_copy['w_median'] = 0
df_copy['w_rrms'] = 0
df_copy['D_mean'] = 0
df_copy['D_std'] = 0
df_copy['D_median'] = 0
df_copy['D_rrms'] = 0

DEFAULT_RETURN = (pd.NA, pd.NA, pd.NA, pd.NA, pd.NA, pd.NA, pd.NA, pd.NA)

GAIA_ID_COLNAME = 'source_id'
SDSS_ID_COLNAME = 'specObjId'

MIN_WVL = 400  # nm
MAX_WVL = 900  # nm


Q_MEDIAN = np.load(os.path.join(DATA_PATH, 'q_median.npy'))


def get_w(*, gaia_id=None, sdss_id=None, k=0.5):
    """
    Gets w for a given SDSS ID.

    """
    if gaia_id is None and sdss_id is None:
        print('No ID supplied')
        return DEFAULT_RETURN
    elif sdss_id is None:
        try:
            sdss_id = DF.loc[DF[GAIA_ID_COLNAME] == gaia_id][SDSS_ID_COLNAME].to_list()[0]
        except:
            print('Could not find an SDSS ID corresponding to the provided Gaia ID')
            return DEFAULT_RETURN
    elif gaia_id is None:
        try:
            gaia_id = DF.loc[DF[SDSS_ID_COLNAME] == sdss_id][GAIA_ID_COLNAME].to_list()[0]
        except:
            print('Could not find an Gaia ID corresponding to the provided SDSS ID')
            return DEFAULT_RETURN

    # obtain the q to correct it for the median error
    q_path = os.path.join(SPECTRA_PATH, f'spectra_G{gaia_id}_S{sdss_id}.npy')

    if not os.path.exists(q_path):
        return DEFAULT_RETURN

    data = decode_npy(q_path)

    sampling = data['sampling']
    mask = (sampling > MIN_WVL) & (sampling < MAX_WVL)

    #print(sampling)
    #print(sampling[mask])
    #print(len(sampling), len(sampling[mask]), len(data['q']))

    sampling = sampling[mask]
    q = data['q']
        
    floors = np.floor(sampling).astype(np.int64)
    closer_to_higher = (sampling - floors) > 0.5
    indices = floors - MIN_WVL + closer_to_higher

    w = q / Q_MEDIAN[indices]

    w_25 = np.quantile(w, 0.25)
    w_75 = np.quantile(w, 0.75)
    w_rrms = 0.7413 * (w_75 - w_25)

    D = np.abs(w - 1)
    D_25 = np.quantile(D, 0.25)
    D_75 = np.quantile(D, 0.75)
    D_rrms = 0.7413 * (D_75 - D_25)

    if not os.path.exists(W_PATH):
        os.makedirs(W_PATH)

    w_path = os.path.join(W_PATH, f'w_G{gaia_id}_S{sdss_id}.npy')
    np.save(w_path, np.array([sampling, w]))

    return (np.mean(w), np.std(w), np.median(w), w_rrms, np.mean(D), np.std(D), np.median(D), D_rrms)


if __name__ == '__main__':
    warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)
    for gaia_id in tqdm(DF[:][GAIA_ID_COLNAME]):
        w_mean, w_std, w_median, w_rrms, D_mean, D_std, D_median, D_rrms = get_w(gaia_id=gaia_id)

        df_copy.loc[df_copy[GAIA_ID_COLNAME] == gaia_id, 'w_mean'] = w_mean
        df_copy.loc[df_copy[GAIA_ID_COLNAME] == gaia_id, 'w_std'] = w_std
        df_copy.loc[df_copy[GAIA_ID_COLNAME] == gaia_id, 'w_median'] = w_median
        df_copy.loc[df_copy[GAIA_ID_COLNAME] == gaia_id, 'w_rrms'] = w_rrms

        df_copy.loc[df_copy[GAIA_ID_COLNAME] == gaia_id, 'D_mean'] = D_mean
        df_copy.loc[df_copy[GAIA_ID_COLNAME] == gaia_id, 'D_std'] = D_std
        df_copy.loc[df_copy[GAIA_ID_COLNAME] == gaia_id, 'D_median'] = D_median
        df_copy.loc[df_copy[GAIA_ID_COLNAME] == gaia_id, 'D_rrms'] = D_rrms

    print('Data processed. Saving to main_table.csv...')
    df_copy.to_csv(os.path.join(DATA_PATH, 'main_table.csv'), index=False)
    print('Done.')
