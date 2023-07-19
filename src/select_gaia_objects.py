import os

from astroquery.gaia import Gaia

from util import get_script_path


query = """
SELECT gs.source_id, sdss.clean_sdssdr13_oid, sdss.original_ext_source_id, sdss.angular_distance,
       gs.ra AS gaia_ra, gs.dec AS gaia_dec, gs.phot_g_mean_mag, gs.bp_g, gs.g_rp
FROM gaiadr3.xp_summary xp
JOIN gaiadr3.sdssdr13_best_neighbour sdss ON sdss.source_id = xp.source_id
JOIN gaiadr3.gaia_source gs ON gs.source_id = xp.source_id
WHERE xp.bp_n_contaminated_transits = 0 AND xp.rp_n_contaminated_transits = 0
AND gs.phot_g_mean_mag < 20.5
AND (gs.ra < 60 OR gs.ra > 300) AND dec > -1.67 AND dec < 1.67
"""

data_path = os.path.join(get_script_path(), '..', 'data')

if not os.path.exists(data_path):
    os.mkdir(data_path)
    print('Created data directory.')

print('Querying Gaia database...')

csv_path = os.path.join(data_path, 'gaia_objects.csv')
Gaia.launch_job_async(query).get_results().to_pandas().to_csv(csv_path, index=False)

print('Done.')
