#!/usr/bin/env python
#-*- coding: UTF-8 -*-

"""
This is a simple implementation of the Fourier-Mellin Trasform
image = [m, n]
ft_magnitude = |fft(image)|
lp_ft_magnitude = logpolar(ft_magnitude)
fmt = fft(lp_ft_magnitude)
"""

from automask import get_mask
from autopipe import showimage
from cache import Cache
from enhance import get_intensity, get_centered, logscale, equalize
from numpy import fft
from numpy import ndarray as atype
from numpy import sin, cos, exp, log, arctan2
from scipy import misc, signal
from scipy.ndimage import geometric_transform
import cv2.cv as cv
import numpy as np

tau = 2 * np.pi


#@Cache("fmt.correlate2d.pickle")
def correlate2d(array1, array2):
    """
    Performs cross correlation between array1 and array2.
    Returns a array of shape h1*2+h2-2, w1*2+w2-2
    """
    rows1, cols1 = array1.shape
    rows2, cols2 = array2.shape
    minrows = min(rows1, rows2)
    mincols = min(cols1, cols2)
    marginrows = (max(rows1, rows2) - 1) / 2.
    margincols = (max(cols1, cols2) - 1) / 2.
    matrix1 = cv.fromarray(np.float32(array1))

    correlation_shape = rows1 * 2 + rows2 - 2, cols1 * 2 + cols2 - 2
    correlation = np.zeros(correlation_shape)
    correlation[rows1 - 1:rows1 - 1 + rows2,cols1-1:cols1-1+cols2] = array2
    correlation_matrix = cv.fromarray(np.float32(correlation))

    result = np.zeros((rows1 + rows2 - 1, cols1 + cols2 - 1))
    result_matrix = cv.fromarray(np.float32(result))

    cv.MatchTemplate(matrix1, correlation_matrix, result_matrix,
        cv.CV_TM_CCORR_NORMED)
    result = np.asarray(result_matrix)
    result = result[marginrows:marginrows + minrows, margincols:margincols + mincols]
    return result


@Cache("fmt.ifft.pickle")
def get_shiftedifft(array):
    shiftedifft = fft.ifftshift(fft.ifft2(array))
    return shiftedifft


@Cache("fmt.fft.pickle")
def get_shiftedfft(array):
    shiftedfft = fft.fftshift(fft.fft2(array))
    return shiftedfft


@Cache("fmt.logpolar.pickle")
def get_logpolar(array, interpolation=0, reverse=False):
    """
    Returns a new array with the logpolar transfamation of array.
    Interpolation can be:
        0 Near
        1 Linear
        2 Bilineal
        3 Cubic
        4
        5
    """
    assert interpolation in range(6)
    rows, cols = array.shape
    row0 = rows / 2.
    col0 = cols / 2.
    theta_scalar = tau / cols
    max_radius = (row0 ** 2 + col0 ** 2) ** .5
    rho_scalar = log(max_radius) / cols

    def cart2logpol(dst_coords):
        theta, rho = dst_coords
        rho = exp(rho * rho_scalar)
        theta = np.pi / 2 - theta * theta_scalar
        row_from = rho * cos(theta) + row0
        col_from = rho * sin(theta) + col0
        return row_from, col_from

    def logpol2cart(dst_coords):
        xindex, yindex = dst_coords
        x = xindex - col0
        y = yindex - row0

        r = np.log(np.sqrt(x ** 2 + y ** 2)) / rho_scalar
        theta = np.arctan2(y, x)
        theta_index = np.round((theta + np.pi) * cols / tau)
        return theta_index, r

    trans = logpol2cart if reverse else cart2logpol

    logpolar = geometric_transform(array, trans, array.shape,
        order=interpolation)
    return logpolar



@Cache("fmt.hi_pass_filter.pickle")
def hi_pass_filter(array, radius=0.2, softness=4):
    radius = round(min(array.shape) * radius)
    window = np.kaiser(radius, softness)
    mask = np.ones_like(array) - get_mask(array.shape, window)
    masked = array * mask
    return masked
    

@Cache("fmt.get_fmt.pickle")
def get_fmt(array):
    """
    Follows this algoritm:
        * FFT with centered frecuencies
        * convolucionar la magnitud con un filtro high pass
        * Logpolar
        * FFT with centered frecuencies
    """
    fourier = get_shiftedfft(array)
    magnitude = np.abs(fourier)
    hi_passed = hi_pass_filter(magnitude, .15, 2)
    logpolar = get_logpolar(hi_passed, 3)
    fmt = get_shiftedfft(logpolar)
    return fmt


def get_correlation(image1, image2):
    """
    Todo esto es un invento y debe ser revisado
    """
    return get_fmt_correlation(image1, image2)
    min_rows = min(image1.shape[0], image2.shape[0])
    min_cols = min(image1.shape[1], image2.shape[1])
    image1 = misc.imresize(image1, (min_rows, min_cols))
    image2 = misc.imresize(image2, (min_rows, min_cols))
    fmt1 = get_fmt(image1)
    fmt2 = get_fmt(image2)
    intensity1 = get_intensity(fmt1)
    intensity2 = get_intensity(fmt2)
    intensitydiff = (intensity2 - intensity1) ** 2
    diff = intensitydiff.mean()
    correlation = (54**2 / (1 + diff)) ** 2
    return correlation


@Cache("fmt.get_fmt_correlation.pickle", 60)
def get_fmt_correlation(image1, image2):
    image1 = get_centered(image1)
    image2 = get_centered(image2)

#    min_rows = min(image1.shape[0], image2.shape[0])
#    min_cols = min(image1.shape[1], image2.shape[1])
#    image1 = misc.imresize(image1, (min_rows, min_cols))
#    image2 = misc.imresize(image2, (min_rows, min_cols))

    fmt1 = get_fmt(image1)
    fmt2 = get_fmt(image2)
    correlation = correlate2d(get_intensity(fmt1), get_intensity(fmt2))
    showimage(equalize(get_shiftedifft(correlation)))
    argmax = np.unravel_index(correlation.argmax(), correlation.shape)
    peak = correlation[argmax]
    relrow = argmax[0] - correlation.shape[0] / 2.
    relcol = argmax[1] - correlation.shape[1] / 2.
    return peak, (relrow, relcol)
