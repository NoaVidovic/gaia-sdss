import os 

from SciServer import Authentication, CasJobs

from util import get_script_path


Authentication.login()

query = """
SELECT objId, specObjId, z, zErr, ra AS specRa, dec AS specDec, photoRa, photoDec,
       psfMag_u, psfMag_g, psfMag_r, psfMag_i, psfMag_z,
       psfMagErr_u, psfMagErr_g, psfMagErr_r, psfMagErr_i, psfMagErr_z
FROM SpecPhoto
WHERE psfMag_r < 20
AND ra > 60 AND ra < 300 AND dec > -1.67 AND dec < 1.67
"""

data_path = os.path.join(get_script_path(), '..', 'data')

if not os.path.exists(data_path):
    os.mkdir(data_path)
    print('Created data directory.')

print("Querying SDSS database...")

csv_path = os.path.join(data_path, 'sdss_objects.csv')
CasJobs.executeQuery(query, context="DR18").to_csv(csv_path, index=False) 

print("Done.")
