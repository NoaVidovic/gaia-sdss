import pandas as pd
import sdss

from astroquery.sdss import SDSS as aq_sdss

from util import get_script_path


DF = pd.read_csv(f'{get_script_path()}/../data/gaia_sdss_joined.csv')

for sdss_id in DF[:1000]['specObjId']:
    try:
        sp = sdss.SpecObj(int(sdss_id))
        data = aq_sdss.get_spectra(plate=sp.plate, mjd=sp.mjd, fiberID=sp.fiberID)
    except:
        continue

    sdss_sampling = 10 ** data[0][1].data['loglam'] / 10

    with open(f'{get_script_path()}/../data/sdss_sampling.csv', 'a') as f:
        f.write(','.join([str(sdss_id)] + [str(x) for x in sdss_sampling]) + '\n')
