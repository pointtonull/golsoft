#!/usr/bin/env python
#-*- coding: UTF-8 -*-

"""
Simple implementation of the Angular Spectrum Method to reconstruct lensless
holograms
"""

from numpy import exp, sqrt
import numpy as np

tau = 6.283185307179586


def get_propagation_array(shape, distance, wavelength, (dx, dy)):
    """
    dx an dy are the phisical distance between pixels on the sensor surface
    """
    wavenumber = tau / wavelength
    rows, cols = shape
    maxrow = rows / 2.
    maxcol = cols / 2.
    minrow, mincol = -maxrow, -maxcol
    row, col = np.ogrid[minrow:maxrow:1, mincol:maxcol:1]
    frow = 1. / (rows * dx) 
    fcol = 1. / (cols * dy)
    phase_correction_factor = wavenumber * sqrt(1 - (wavelength * frow *
        row) ** 2 - (wavelength * fcol * col) ** 2)
    propagation_array = exp(1j * phase_correction_factor * distance)
    return propagation_array
