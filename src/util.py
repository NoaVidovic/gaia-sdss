import numpy as np
import os
import sys


def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))


def decode_npy(npy_path):
    sampling, sdss_flux, sdss_flux_fit, sdss_conv, gaia_flux, q = np.load(npy_path, allow_pickle=True)
    return {'sampling': sampling, 'sdss_flux': sdss_flux, 'sdss_flux_fit': sdss_flux_fit, 'sdss_conv': sdss_conv,
            'gaia_flux': gaia_flux, 'q': q}

def decode_npy_w(npy_path):
    sampling, w = np.load(npy_path, allow_pickle=True)
    return {'sampling': sampling, 'w': w}
