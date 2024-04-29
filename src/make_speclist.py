import sdss
import pandas as pd


def specid_to_name(spec_id):
    try:
        sp = sdss.SpecObj(spec_id)
        return f'{sp.plate:04}/spec-{sp.plate:04}-{sp.mjd:05}-{sp.fiberID:04}.fits'
    except:
        return None


###

spec_obj_ids = pd.read_csv('data/main_from_gaia_archive.csv')['specObjId']
speclist = spec_obj_ids.map(specid_to_name).dropna()
speclist.to_csv('data/speclist.csv')
