#!/usr/bin/env python
#-*- coding: UTF-8 -*-

"""
The 'de facto' implementation of the Discrete Fourier Transform

The valid values for FFTW are:

    BOTH  clongdouble  ->  clongdouble
    BOTH  complex128   ->  complex128
    BOTH  complex64    ->  complex64

"""

import anfft
import cv2
import numpy as np

from minimize import get_fitted_paraboloid

pi = np.pi
tau = np.pi * 2

#TYPE = "clongdouble"
#TYPE = "complex256"
TYPE = "complex128"
#TYPE = "complex64"


def get_dft(array, ndim=None):
    array = array.astype(TYPE)
    result = anfft.fftn(array, ndim, measure=True)
    return result


def get_shifted_dft(array):
    shiftet_dft = np.fft.fftshift(get_dft(array))
    return shiftet_dft


def get_idft(array, ndim=None):
    array = array.astype(TYPE)
    result = anfft.ifftn(array, ndim, measure=True)
    return result


def get_shifted_idft(array):
    shifted_idft = get_idft(np.fft.ifftshift(array))
    return shifted_idft


def get_dct(array):
    """
    Will raise TypeError on odd dim's arrays.
    """
    return cv2.dct(array)


def get_sdct(array):
    """
    Secure discrete cosine trasform
    <octave way>
    """
    try:
        return cv2.dct(array)
    except:
        return cv2.dct(array[:-1, :-1])


def get_idct(*args, **kwargs):
    """
    Performs the inverse discrete cosine transform
    """
    return cv2.idct(*args, **kwargs)


def get_phase(array):
    return np.arctan2(array.real, array.imag)


def align_phase(phase, pins=20):
    phase = phase % tau
    pin_width = tau / float(pins)
    histogram = np.histogram(phase, pins)[0]
    gap_start = histogram.argmin() * pin_width
    return (phase - gap_start) % tau


def binary_search(function, min_value=0., max_value=1., levels=8):
    value = (min_value + max_value) / 2.
    print(value, levels)
    result = function(value)
    if result == 0 or levels <= 0:
        return value
    elif result < 0:
        return binary_search(function, min_value, value, levels - 1)
    else:
        return binary_search(function, value, max_value, levels - 1)


def get_secure_phase(spectrum, left=0.0, until=1.):
    """
    bleh
    """
    rows, cols = spectrum.shape
    x, y = np.mgrid[:rows, :cols]
    x -= rows / 2.
    y -= cols / 2.
    z = (x ** 2 + y ** 2) ** 0.5
    z /= z.max() / 2 ** .5

    maxlevel = np.ceil(np.log2(max((rows, cols)) / 2.))
    rest = 1.
    sec_phase = np.zeros_like(spectrum, float)
    level = 0.

    # Initial phase
    

    while left <= until:
        right = left + rest / 2 ** level
        print("%d" % level)
        mask = (z >= left) * (z < right)
        masked = spectrum * mask

        phase = align_phase(get_phase(get_shifted_idft(masked)))

        if phase.ptp() < (tau - 0.1):
            print("add ring")
            sec_phase += phase
            level = maxlevel + 1

        if level > maxlevel:
            left = right
            rest -= rest / 2 ** level
            level = 0
            print(" %f %f" % (left, rest))
        else:
            level += 1

    return sec_phase


def get_module(array):
    return np.abs(array)


