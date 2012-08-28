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
LAMBDA = 6.328e-07 # wave length
DX = 8.39e-6
DY = 8.46e-6
K = tau / LAMBDA # wave number
EPSILON = 1e-16


@cache.hybrid(reset=0)
def get_propagation_array(shape, distance):
    rows, cols = shape
    maxrow = rows / 2.
    maxcol = cols / 2.
    minrow, mincol = -maxrow, -maxcol
    row, col = np.ogrid[minrow:maxrow:1, mincol:maxcol:1]
    frow = 1. / (rows * DX) 
    fcol = 1. / (cols * DY)
    phase_correction_factor = K * sqrt(1 - (LAMBDA * frow * row) ** 2
        - (LAMBDA * fcol * col) ** 2)
    propagation_array = exp(1j * phase_correction_factor * distance)
    return propagation_array


