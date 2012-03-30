#!/usr/bin/env python
#-*- coding: UTF-8 -*-

"""
This is a simple implementation of the Fourier-Mellin Trasform
"""

from automask import get_mask
from autopipe import showimage
from image import get_intensity, get_centered, logscale, equalize
from numpy import fft
from numpy import ndarray as atype
from numpy import sin, cos, exp, log, arctan2
from scipy import misc, signal
from scipy.ndimage import geometric_transform
import cache
import cv2.cv as cv
import numpy as np

tau = 2 * np.pi


@cache.hybrid
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
    result = result[marginrows:marginrows + minrows,
        margincols:margincols + mincols]
    return result


@cache.hybrid
def get_fft(array):
    fft = np.fft.fft2(array)

#    fft = np.zeros(array.shape)
#    fft_matrix = cv.fromarray(fft)
#    matrix = cv.fromarray(np.float32(array))
#    cv.DFT(matrix, fft_matrix, 0)
#    fft = np.asarray(fft_matrix)
#    showimage(equalize(fft))

    return fft


def get_shiftedfft(array):
    shiftedfft = fft.fftshift(get_fft(array))
    return shiftedfft


@cache.hybrid
def get_ifft(array):
    return np.fft.ifft2(array)


def get_shiftedifft(array):
    shiftedifft = fft.ifftshift(get_ifft(array))
    return shiftedifft


@cache.hybrid
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



@cache.hybrid
def hi_pass_filter(array, radius=0.2, softness=4):
    radius = round(min(array.shape) * radius)
    window = np.kaiser(radius, softness)
    mask = np.ones_like(array) - get_mask(array.shape, window)
    masked = array * mask
    return masked
    

@cache.hybrid
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


@cache.hybrid
def get_fmt_correlation(image1, image2):
    image1 = get_centered(image1)
    image2 = get_centered(image2)

    #1 Apply same 2D windowing to both the images
    #2 Fourier-Mellin transform (which is translation invariant:
    #  is a 2D Fourier transform followed by abs()) (sic)
    #3 map the transformed images to log-polar space, so that rotation / scaling
    #  become Dx/Dy translations on the new axes
    #4 take 2D Fourier transform to results
    fmt1 = get_fmt(image1)
    fmt2 = get_fmt(image2)

    # Multiply one for the coniugate of the other
    #5 divide (element by element) this product by its abs()

#    correlation = correlate2d(get_intensity(fmt1), get_intensity(fmt2))
    correlation = correlate2d(np.abs(fmt1), np.abs(fmt2)) #magnitude
#    correlation = correlate2d(np.angle(fmt1), np.angle(fmt2)) #phase
    showimage(logscale(correlation), equalize(correlation))

    #6 reverse 2D Fourier transform of result of step 5; note that steps 4+6 are
    #  equivalent to perform a correlation on results of step 3; step 5
    #  introduces a normalization; all together is the "cross power spectrum"
    #  calculation

    #7 calculate the position of peaks; after appropriate conversion these are
    #  candidates for rotation and scale parameters. Note that:
    #  a. for scale it will be probably necessary to test more than one peak
    #  b. rotation is affected by a 180 degree ambiguity therefore at least two 
    #     cases shell be tested.
    #8 adjust images for scale and rotation
    #9 find Dx/Dy and adjust for translation (more or less the same steps 4 to 7,
    #  but applied to results of step 8 and with less problems).

    argmax = np.unravel_index(correlation.argmax(), correlation.shape)
    peak = correlation[argmax]
    relrow = argmax[0] - correlation.shape[0] / 2.
    relcol = argmax[1] - correlation.shape[1] / 2.
    return peak, (relrow, relcol)
