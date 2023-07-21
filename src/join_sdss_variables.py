import os
import math
import pandas as pd

from util import get_script_path

data_path = os.path.join(get_script_path(), '..', 'data')

gaia_sdss_joined = pd.read_csv(os.path.join(data_path, 'gaia_sdss_joined.csv'))
variables = pd.read_csv(os.path.join(data_path, 'variables.csv'))

def find_nearby_stars(df1, df2):
    i=1
    row_count = df1.shape[0]
    # Loop through stars in file1 and file2 to calculate distances
    for _, star1 in df1.iterrows():
        for _, star2 in df2.iterrows():

            distance = haversine_distance(star1['photoRa'], star1['photoDec'], star2['ra'], star2['dec'])

            if distance < 1:

                df1.loc[df1.index == star1.name, 'variablesID'] = int(star2['ID'])
                df1['variablesID'] = pd.to_numeric(df1['variablesID'], errors='coerce').astype(pd.Int64Dtype())
                df1.to_csv(os.path.join(data_path, 'gaia_sdss_joined_with_variablesID.csv'), index=False)
                print('found one')
                continue

        percentage = (i/row_count) * 100
        print(str(round(percentage, 2)) + '% Done')
        i += 1
    return


def haversine_distance(ra1, dec1, ra2, dec2):
    # Convert degrees to radians
    lat_rad1 = dec1 * (math.pi/180)
    lon_rad1 = ra1 * (math.pi/180)

    lat_rad2 = dec2 * (math.pi/180)
    lon_rad2 = ra2 * (math.pi/180)

    # Calculate the Haversine distance in radians
    d_rad = 2 * math.asin(math.sqrt(math.sin((lat_rad2 - lat_rad1)/2)**2 +
                                    math.cos(lat_rad1) * math.cos(lat_rad2) * math.sin((lon_rad2 - lon_rad1)/2)**2))

    # Convert radians to degrees
    d_deg = d_rad * (180 / math.pi)

    # Convert degrees to arcseconds
    d_arcsec = d_deg * 3600

    return d_arcsec


find_nearby_stars(gaia_sdss_joined, variables)