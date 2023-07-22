import os
import pandas as pd
import numpy as np
import sys
from sklearn.neighbors import BallTree

def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))


data_path = os.path.join(get_script_path(), '..', 'data')

gaia_sdss_joined = pd.read_csv(os.path.join(data_path, 'gaia_sdss_joined.csv'))
variables = pd.read_csv(os.path.join(data_path, 'variables.csv'))

# Convert degrees to radians for calculations
gaia_sdss_joined[['photoRa', 'photoDec']] = np.radians(gaia_sdss_joined[['photoRa', 'photoDec']])
variables[['ra', 'dec']] = np.radians(variables[['ra', 'dec']])


def find_nearby_stars(df1, df2):
    # Create a BallTree for efficient nearest neighbor search
    tree = BallTree(df2[['ra', 'dec']], metric='haversine')

    # Query the tree for neighbors within given radius
    dist = tree.query_radius(df1[['photoRa', 'photoDec']], r=0.1/3600)

    df1['variablesID'] = pd.NA
    variable_multiples = pd.DataFrame(columns=['source_id', 'variablesID'])

    for i, matches in enumerate(dist):
        if len(matches) == 1:  # One match
            df1.loc[i, 'variablesID'] = df2.loc[matches[0], 'ID']
        elif len(matches) > 1:  # Multiple matches
            new_rows = [{'source_id': df1.loc[i, 'source_id'], 'variablesID': df2.loc[match, 'ID']} for match in matches]
            variable_multiples = pd.concat([variable_multiples, pd.DataFrame(new_rows)], ignore_index=True)

    df1.to_csv(os.path.join(data_path, 'gaia_sdss_joined_with_variablesID.csv'), index=False)
    variable_multiples.to_csv(os.path.join(data_path, 'variable_multiples.csv'), index=False)

    return df1, variable_multiples


find_nearby_stars(gaia_sdss_joined, variables)
