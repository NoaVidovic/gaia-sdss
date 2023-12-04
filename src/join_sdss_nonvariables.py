import os
import pandas as pd
import numpy as np
import sys
from sklearn.neighbors import BallTree


def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))


data_path = os.path.join(get_script_path(), '..', 'data')

gaia_sdss_joined = pd.read_csv(os.path.join(data_path, 'gaia_sdss_joined_with_variablesID.csv'))
nonvariables = pd.read_csv(os.path.join(data_path, 'nonvariables.csv'))

# Convert degrees to radians for calculations
nonvariables[['ra', 'dec']] = np.radians(nonvariables[['RA', 'Dec']])


def find_nearby_stars(df1, df2):
    # Create a BallTree for efficient nearest neighbor search
    tree = BallTree(df2[['ra', 'dec']], metric='haversine')

    # Query the tree for neighbors within given radius
    dist = tree.query_radius(df1[['photoRa', 'photoDec']], r=np.radians(0.1/3600))

    df1['nonvariablesID'] = pd.NA
    nonvariable_multiples = pd.DataFrame(columns=['source_id', 'nonvariablesID'])

    for i, matches in enumerate(dist):
        if len(matches) == 1:  # One match
            df1.loc[i, 'nonvariablesID'] = df2.loc[matches[0], 'ID']
        elif len(matches) > 1:  # Multiple matches
            new_rows = [{'source_id': df1.loc[i, 'source_id'], 'nonvariablesID': df2.loc[match, 'ID']} for match in matches]
            nonvariable_multiples = pd.concat([nonvariable_multiples, pd.DataFrame(new_rows)], ignore_index=True)

    df1.to_csv(os.path.join(data_path, 'gaia_sdss_joined_nonvariablesID.csv'), index=False)
    nonvariable_multiples.to_csv(os.path.join(data_path, 'nonvariable_multiples.csv'), index=False)

    return df1, nonvariable_multiples


find_nearby_stars(gaia_sdss_joined, nonvariables)
