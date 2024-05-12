import os
import pandas as pd
import numpy as np

from astropy.io import fits
from sklearn.neighbors import BallTree
from time import time

from util import get_script_path


ARCSECOND = 1/3600 * np.pi/180  # in radians


def dat_line_to_entry(l):
    _id, ra, dec = l.split()[:3]
    ra = float(ra) * np.pi / 180
    dec = float(dec) * np.pi / 180

    return [ _id, ra, dec ]


def read_dat(path):
    with open(path, 'r') as f:
        entries = [ dat_line_to_entry(l) for l in f.readlines() if l[0] != '#' ]

    return pd.DataFrame(entries, columns=['id', 'ra', 'dec'])


def read_fits(path):
    qf = fits.open(path)[1].data

    q_id = qf.field('OBJID').astype(str)
    q_ra = qf.field('RA').astype(np.float64) * np.pi/180
    q_dec = qf.field('DEC').astype(np.float64) * np.pi/180

    return pd.DataFrame(np.array([q_id, q_ra, q_dec]).transpose(), columns=['id', 'ra', 'dec'])


def xmatch(main_df, ext_filename):
    name, ext = ext_filename.split('.')
    if ext == 'fits':
        ext_df = read_fits(os.path.join(data_path, ext_filename))
    elif ext == 'dat':
        ext_df = read_dat(os.path.join(data_path, ext_filename))
    else:
        print('Bad extension')
        return None

    # Create a BallTree for efficient nearest neighbor search
    tree = BallTree(ext_df[['dec', 'ra']], metric='haversine')

    main_rad = np.array((main_df.photoDec * np.pi/180, main_df.photoRa * np.pi/180)).transpose()

    # Query the tree for neighbors within 0.1 as
    neighbor_list = tree.query_radius(main_rad, r=0.1*ARCSECOND)

    main_df[f'{name}_match_id'] = 'NA'
    main_df[f'{name}_match_num'] = pd.NA

    for i, neighbors in enumerate(neighbor_list):
        main_df.loc[i, f'{name}_match_num'] = len(neighbors)
        if len(neighbors) == 1:
            main_df.loc[i, f'{name}_match_id'] = ext_df.loc[neighbors[0], 'id']

    return main_df


data_path = os.path.join(get_script_path(), '..', 'data')
main_table = pd.read_csv(os.path.join(data_path, 'main_table_S82.csv'))

print('Cross-matching with non-variable star dataset...')
t0 = time()
main_table = xmatch(main_table, 'nonvar.dat')
print(f'Done in {time()-t0:.3f} seconds.')

print('Cross-matching with variable star dataset...')
t0 = time()
main_table = xmatch(main_table, 'var.dat')
print(f'Done in {time()-t0:.3f} seconds.')

print('Cross-matching with quasar dataset...')
t0 = time()
main_table = xmatch(main_table, 'quasar.fits')
print(f'Done in {time()-t0:.3f} seconds.')

print('Saving to file...')
main_table.to_csv(os.path.join(data_path, 'main_table_S82.csv'), index=False)
print('Done.')
