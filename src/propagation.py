#!/usr/bin/env python
#-*- coding: UTF-8 -*-

"""
Simple implementation of the Angular Spectrum Method to reconstruct lensless
holograms
"""

from numpy import exp, sqrt
import numpy as np

import cache


tau = 6.283185307179586
DX = 8.39e-6
DY = 8.46e-6
EPSILON = 1e-16


@cache.hybrid(reset=0)
def get_propagation_array(shape, distance, wavelength):
    wavenumber = tau / wavelength
    rows, cols = shape
    maxrow = rows / 2.
    maxcol = cols / 2.
    minrow, mincol = -maxrow, -maxcol
    row, col = np.ogrid[minrow:maxrow:1, mincol:maxcol:1]
    frow = 1. / (rows * DX) 
    fcol = 1. / (cols * DY)
    phase_correction_factor = wavenumber * sqrt(1 - (wavelength * frow *
        row) ** 2 - (wavelength * fcol * col) ** 2)
    propagation_array = exp(1j * phase_correction_factor * distance)
    return propagation_array
