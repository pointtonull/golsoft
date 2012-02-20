#!/usr/bin/env python
#-*- coding: UTF-8 -*-

"""
This is a simple implementation of the Fourier-Mellin Trasform
image = [m, n]
ft_magnitude = |fft(image)|
lp_ft_magnitude = logpolar(ft_magnitude)
fmt = fft(lp_ft_magnitude)
"""

from numpy import pi, sin, cos, exp, ndarray, asarray, log
from scipy.ndimage import geometric_transform

tau = 2 * pi

def get_logpolar(array):
    rows, cols = array.shape
    out_rows = rows
    out_cols = cols
    row0 = rows / 2.
    col0 = cols / 2.
    theta_scalar = tau / out_cols
    max_radius = (row0 ** 2 + col0 ** 2) ** .5
    rho_scalar = log(max_radius) / cols

    def out2in(dst_coords):
        theta, rho = dst_coords
        rho = exp(rho * rho_scalar)
        theta = theta * theta_scalar - pi
        row_from = rho * cos(theta) + row0
        col_from = rho * sin(theta) + col0
        return row_from, col_from

    logpolar = geometric_transform(array, out2in, (out_rows, out_cols), order=3)
    return logpolar
