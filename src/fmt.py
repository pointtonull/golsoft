#!/usr/bin/env python
#-*- coding: UTF-8 -*-

"""
This is a simple implementation of the Fourier-Mellin Trasform
"""

from autopipe import showimage
from automask import get_mask
from dft import get_shifted_dft
from image import get_centered, get_logpolar
import cache
import cv2.cv as cv
import numpy as np

tau = 6.283185307179586


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
    correlation[rows1 - 1:rows1 - 1 + rows2, cols1-1:cols1 - 1 + cols2] = array2
    correlation_matrix = cv.fromarray(np.float32(correlation))

    result = np.zeros((rows1 + rows2 - 1, cols1 + cols2 - 1))
    result_matrix = cv.fromarray(np.float32(result))

    cv.MatchTemplate(matrix1, correlation_matrix, result_matrix,
        cv.CV_TM_CCORR_NORMED)
    result = np.asarray(result_matrix)
    result = result[marginrows:marginrows + minrows,
        margincols:margincols + mincols]
    return result


@cache.hybrid(reset=0)
def get_high_pass_filter(array, radius=0.2, softness=4):
    radius = round(min(array.shape) * radius)
    window = np.kaiser(radius, softness)
    mask = np.ones_like(array) - get_mask(array.shape, window)
    masked = array * mask
    return masked
    

@cache.hybrid
def get_fmt(array, high_pass=0.15):
    """
    Follows this algoritm:
        * FFT with centered frecuencies
        * convolucionar la magnitud con un filtro high pass
        * Logpolar
        * FFT with centered frecuencies
    """
    fourier = get_shifted_dft(array)
    magnitude = np.abs(fourier)
    high_passed = get_high_pass_filter(magnitude, high_pass, 0)
    logpolar = get_logpolar(high_passed, 3)
    fmt = get_shifted_dft(logpolar)
    return fmt


@cache.hybrid
def get_fmt_correlation(image1, image2, high_pass=0.15):
    image1 = get_centered(image1)
    image2 = get_centered(image2)

    fmt1 = get_fmt(image1, high_pass)
    fmt2 = get_fmt(image2, high_pass)

    correlation = correlate2d(np.abs(fmt1), np.abs(fmt2)) #magnitude

    argmax = np.unravel_index(correlation.argmax(), correlation.shape)
    peak = correlation[argmax]
    relrow = argmax[0] - correlation.shape[0] / 2.
    relcol = argmax[1] - correlation.shape[1] / 2.
    return peak, (relrow, relcol)
