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
