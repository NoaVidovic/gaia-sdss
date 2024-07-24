import os
import pandas as pd

from util import get_script_path


data_path = os.path.join(get_script_path(), '..', 'data')

print('Joining Gaia standard deviation data...')
main = pd.read_csv(os.path.join(data_path, 'main_table.csv'))
gaia_std = pd.read_csv(os.path.join(data_path, 'gaia_std.csv'))

pd.merge(main, gaia_std, on='source_id', how='left').to_csv(os.path.join(data_path, 'main_table.csv'), index=False)

print('Done.')
