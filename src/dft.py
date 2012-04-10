#!/usr/bin/env python
#-*- coding: UTF-8 -*-

"""
The 'de facto' implementation of Fast Fourier Transform

The valid values for FFTW are:

    BOTH  clongdouble  ->  clongdouble
    BOTH  complex128   ->  complex128
    BOTH  complex64    ->  complex64

"""

from autopipe import showimage
from image import equalize
import anfft
import cache
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
