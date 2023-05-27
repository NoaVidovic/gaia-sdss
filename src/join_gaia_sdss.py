import os
import pandas as pd

from util import get_script_path


data_path = os.path.join(get_script_path(), '..', 'data')

print('Joining Gaia and SDSS data...')
gaia = pd.read_csv(os.path.join(data_path, 'gaia_objects.csv'))
sdss = pd.read_csv(os.path.join(data_path, 'sdss_objects.csv'))

pd.merge(gaia, sdss, left_on='original_ext_source_id', right_on='objId').to_csv(os.path.join(data_path, 'gaia_sdss_joined.csv'), index=False)

print('Done.')
