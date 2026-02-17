# ShapePipe code for compiling PSFs

"""PSFEX INTERPOLATION SCRIPT.

This module computes the PSFs from a PSFEx model at several galaxy positions.

:Authors: Morgan Schmitz and Axel Guinot

"""

import os
import re

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval
#from sqlitedict import SqliteDict

#from shapepipe.pipeline import file_io

try:
    import galsim.hsm as hsm
    from galsim import Image
except ImportError:
    import_fail = True
else:
    import_fail = False


NOT_ENOUGH_STARS = 'Fail_stars'
BAD_CHI2 = 'Fail_chi2'
FILE_NOT_FOUND = 'File_not_found'

def interpsfex(dotpsfpath, pos, thresh_star, thresh_chi2):
    """Interpolate PSFEx.

    Use PSFEx generated model to perform spatial PSF interpolation.

    Parameters
    ----------
    dotpsfpath : str
        Path to ``.psf`` file (PSFEx output)
    pos : numpy.ndarray
        Positions where the PSF model should be evaluated
    thresh_star : int
        Threshold of stars under which the PSF is not interpolated
    thresh_chi2 : int
        Threshold for chi squared

    Returns
    -------
    numpy.ndarray
        Array of PSFs, each row is the PSF image at the corresponding position
        requested

    """
    if not os.path.exists(dotpsfpath):
        return FILE_NOT_FOUND

    # read PSF model and extract basis and polynomial degree and scale position
    PSF_model = fits.open(dotpsfpath)[1]

    # Check number of stars used to compute the PSF
    if PSF_model.header['ACCEPTED'] < thresh_star:
        return NOT_ENOUGH_STARS
    if PSF_model.header['CHI2'] > thresh_chi2:
        return BAD_CHI2

    PSF_basis = np.array(PSF_model.data)[0][0]
    try:
        deg = PSF_model.header['POLDEG1']
    except KeyError:
        # constant PSF model
        return PSF_basis[0, :, :]

    # scale coordinates
    x_interp, x_scale = (
        PSF_model.header['POLZERO1'],
        PSF_model.header['POLSCAL1']
    )
    y_interp, y_scale = (
        PSF_model.header['POLZERO2'],
        PSF_model.header['POLSCAL2']

    )
    xs, ys = (pos[:, 0] - x_interp) / x_scale, (pos[:, 1] - y_interp) / y_scale

    # compute polynomial coefficients
    coeffs = np.array([[x ** idx for idx in range(deg + 1)] for x in xs])
    cross_coeffs = np.array([
        np.concatenate([
            [(x ** idx_j) * (y ** idx_i) for idx_j in range(deg - idx_i + 1)]
            for idx_i in range(1, deg + 1)
        ])
        for x, y in zip(xs, ys)
    ])
    coeffs = np.hstack((coeffs, cross_coeffs))

    # compute interpolated PSF
    PSFs = np.array([
        np.sum(
            [coeff * atom for coeff, atom in zip(coeffs_posi, PSF_basis)],
            axis=0,
        )
        for coeffs_posi in coeffs
    ])

    return PSFs


# PSF locations
f070w_pos = np.array([[711,461]])
f480m_pos = np.array([[2945,2705]])
f070w_psfpath = 'psfcompilation_J0025/PSF_data/flux_f070w_img.fluxpsf.fits'
f480m_psfpath = 'psfcompilation_J0025/PSF_data/flux_f480m_img.fluxpsf.fits'
f070w_result = interpsfex(f070w_psfpath, f070w_pos, 0, 1000)
f480m_result = interpsfex(f480m_psfpath, f480m_pos, 0, 1000)


# save to .FITS file
f070w_hdu = fits.PrimaryHDU(data=f070w_result[0])
f070w_hdul = fits.HDUList([f070w_hdu])
f070w_hdul.writeto('F070Wpsf.fits')

f480m_hdu = fits.PrimaryHDU(data=f480m_result[0])
f480m_hdul = fits.HDUList([f480m_hdu])
f480m_hdul.writeto('F480Mpsf.fits')