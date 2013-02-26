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
from matplotlib.pyplot import plot, show

from minimize import get_fitted_paraboloid
from autopipe import showimage

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
    """ Performs the inverse discrete cosine
    transform """
    return cv2.idct(*args, **kwargs)


def get_phase(array):
    return np.arctan2(array.real, array.imag)


def align_phase(phase, pins=20):
    pin_width = tau / float(pins)
    histogram = np.histogram(phase, pins)[0]
    
#    print(histogram) plot(histogram) show(block=True)

    gap_start = histogram.argmin() * pin_width
    phase = (phase - gap_start) % tau

    histogram = np.histogram(phase, pins)[0]
    histogram[histogram < histogram.max() / 100.] = 0

#    print(histogram) plot(histogram) show(block=True)

    is_wrapped = (histogram[0] * histogram[-1]) > 0
    return is_wrapped, phase


class Bisector:
    def __init__(self, universe, left=None, right=None):
        self.universe = universe
        self.left = left or 0
        self.right = right or len(universe)

    def value(self):
        return self.universe[int((self.left + self.right) / 2.)]

    def to_left(self):
        self.right = int(round((self.left + self.right) / 2.))

    def to_right(self):
        self.left = int(round((self.left + self.right) / 2.))

    def is_closed(self):
        return (self.right - self.left) <= 1

    def reset(self):
        self.left = 0
        self.right = len(self.universe)

    def reset_to_right(self):
        self.left = int((self.left + self.right) / 2.) + 1
        self.right = len(self.universe)

    def get_window(self):
        return self.universe[self.left:self.right]

    def __repr__(self):
        return str(self.get_window())


def get_secure_phase(spectrum, until=1.):
    """ If we knew what it was we were doing it would not be called research."""
    rows, cols = spectrum.shape

    def get_partial_phase(max_value):
        mask = spectrum >= max_value
        masked = spectrum * mask
        phase = get_phase(get_shifted_idft(masked))
        paraboloid = get_fitted_paraboloid(phase)
        phase -= paraboloid
        return paraboloid, phase % tau

    # Initial phase

    sec_phase = np.zeros_like(spectrum, float)
    values = sorted(spectrum.ravel(), reverse=True)
    bisector = Bisector(values)
    while bisector.value() > values.min():
        paraboloid, phase = get_partial_phase(bisector.value())
        is_wrapped, phase = align_phase((phase - sec_phase) % tau, 20)
        print(is_wrapped)

        if is_wrapped:
            bisector.to_left()
        else:
            print("add ring")
            sec_phase += phase
            sec_phase += paraboloid
            showimage(sec_phase)
            bisector.reset_to_right()

    return sec_phase


def get_module(array):
    return np.abs(array)


