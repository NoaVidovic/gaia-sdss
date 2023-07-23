import os
import pandas as pd
import numpy as np
from util import get_script_path
from sklearn.neighbors import BallTree

data_path = os.path.join(get_script_path(), '..', 'data')

gaia_sdss_joined = pd.read_csv(os.path.join(data_path, 'gaia_sdss_joined.csv'))
variables = pd.read_csv(os.path.join(data_path, 'variables.csv'))

# Convert degrees to radians for calculations
variables[['ra', 'dec']] = np.radians(variables[['ra', 'dec']])

def find_nearby_stars(df1, df2):
    # Create a BallTree for efficient nearest neighbor search
    tree = BallTree(df2[['ra', 'dec']], metric='haversine')

    # Query the tree for neighbors within given radius
    ind = tree.query_radius(df1[['photoRa', 'photoDec']], r=np.radians(0.1/3600))

    # Initialize the variablesID column with empty values
    df1['variablesID'] = pd.NA

    # Initialize an empty DataFrame for multiple matches
    variable_multiples = pd.DataFrame(columns=['source_id', 'variablesID'])

    # Handle each case (no matches, one match, multiple matches) separately
    for i, matches in enumerate(ind):
        if len(matches) == 1:  # One match
            df1.loc[i, 'variablesID'] = df2.loc[matches[0], 'ID']
        elif len(matches) > 1:  # Multiple matches
            new_rows = [{'source_id': df1.loc[i, 'source_id'], 'variablesID': df2.loc[match, 'ID']} for match in matches]
            variable_multiples = pd.concat([variable_multiples, pd.DataFrame(new_rows)], ignore_index=True)

    # Write the result to a new CSV file
    df1.to_csv(os.path.join(data_path, 'gaia_sdss_joined_with_variablesID.csv'), index=False)

    # Write multiple matches to a new CSV file
    variable_multiples.to_csv(os.path.join(data_path, 'variable_multiples.csv'), index=False)

    return df1, variable_multiples

# Call the function
result, multiples = find_nearby_stars(gaia_sdss_joined, variables)
