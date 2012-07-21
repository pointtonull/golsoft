# !/usr/bin/env python
# -*- coding: UTF-8 -*-

from numpy import arctan2, exp
import numpy as np

tau = 6.283185307179586
pi = 3.1415926589
LAMBDA = 6.328e-07 # wave length
DX = 8.39e-6
DY = 8.46e-6
K = tau / LAMBDA # wave number
EPSILON = 1e-16



def angle2(array):
    return arctan2(array.real, array.imag)


def get_chirp(shape, distance):
    maxrow = shape[0] / 2.
    maxcol = shape[1] / 2.
    minrow, mincol = -maxrow, -maxcol
    row, col = np.ogrid[minrow:maxrow:1., mincol:maxcol:1.]
    fresnel = exp(-1j * pi/(LAMBDA * distance) * (row ** 2 * DX ** 2
        + col ** 2 * DY ** 2))
    return fresnel


def get_wrapped_formula_array(shape, distance):
    rows, cols = shape    
    maxrow = shape[0] / 2.
    maxcol = shape[1] / 2.
    minrow, mincol = -maxrow, -maxcol
    row, col = np.ogrid[minrow:maxrow:1., mincol:maxcol:1.]
    eta = LAMBDA * distance / (rows * DX)
    nu = LAMBDA * distance / (cols * DY)
    wfa = (1j / (LAMBDA * distance)) * exp(-1j * 2 * distance * pi / LAMBDA) * exp(-1j * pi 
        / (LAMBDA * distance) * (row ** 2 * eta ** 2 + col ** 2 * nu ** 2))
    return wfa
