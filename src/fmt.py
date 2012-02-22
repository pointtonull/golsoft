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
from numpy import fft
from numpy import ndarray as atype
from numpy import sin, cos, exp, log
from scipy.ndimage import geometric_transform
import numpy as np

tau = 2 * np.pi

def get_logpolar(array, order=0):
    rows, cols = array.shape
    row0 = rows / 2.
    col0 = cols / 2.
    theta_scalar = tau / cols
    max_radius = (row0 ** 2 + col0 ** 2) ** .5
    rho_scalar = log(max_radius) / cols

    def out2in(dst_coords):
        theta, rho = dst_coords
        rho = exp(rho * rho_scalar)
        theta = theta * theta_scalar - tau / 2
        row_from = rho * cos(theta) + row0
        col_from = rho * sin(theta) + col0
        return row_from, col_from

    logpolar = geometric_transform(array, out2in, array.shape, order=order)
    return logpolar


def hi_pass_filter(array, radius=0.2):
    radius = min(array.shape) * radius
    window = np.bartlett(radius)
    mask = np.ones_like(array) - get_mask(array.shape, window)
    masked = array * mask
    return masked
    

def get_fmt(array):
    """
    Follows this algoritm:
        * FFT with centered frecuencies
        * convolucionar la magnitud con un filtro high pass #TODO
        * Logpolar
        * FFT with centered frecuencies
    """
    fourier = fft.fftshift(fft.fft2(array))
    real_hi_passed = hi_pass_filter(fourier.real, .4)
    imag_hi_passed = hi_pass_filter(fourier.imag, .4)
    real_logpolar = get_logpolar(real_hi_passed, 3)
    imag_logpolar = get_logpolar(imag_hi_passed, 3)
    logpolar = real_logpolar + 1j * imag_logpolar
    fmt = fft.fftshift(fft.fft2(logpolar))
    return fmt
