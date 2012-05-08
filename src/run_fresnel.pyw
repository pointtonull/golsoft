# !/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Testing case
"""

from autopipe import showimage
from image import equalize, imread, normalize, imwrite, get_intensity
from pea import calculate_director_cosines, get_ref_beam, get_propagation_array
from pea import get_auto_mask
from numpy import abs, arctan2, exp
from dft import get_shifted_dft, get_shifted_idft, get_dft
import numpy as np
from ranges import frange
import sys


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
    maxrow = shape[0] / 2.
    maxcol = shape[1] / 2.
    minrow, mincol = -maxrow, -maxcol
    row, col = np.ogrid[minrow:maxrow:1., mincol:maxcol:1.]
    eta = LAMBDA * distance / (512 * DX)
    nu = LAMBDA * distance / (512 * DY)
    wfa = (1j / (LAMBDA * distance)) * exp(-1j * 2 * distance * pi / LAMBDA) * exp(-1j * pi 
        / (LAMBDA * distance) * (row ** 2 * eta ** 2 + col ** 2 * nu ** 2))
    return wfa



def main():
    """
    Le main rutine
    """
    images = [(filename, imread(filename, True)) for filename in sys.argv[1:]]
    for filename, hologram in images:
        print(filename)
        shape = hologram.shape
        hologram = hologram - hologram.mean()
        showimage(equalize(hologram))
        
        distance = .11
        print("Distance: %3.2f" % distance)

        cos_alpha, cos_beta = calculate_director_cosines(hologram)
        print("Cosines: %4.3f %4.3f" % (cos_alpha, cos_beta))

        print("Reference beam: normalized(imag == real)")
        ref_beam = get_ref_beam(shape, cos_alpha, cos_beta)
        showimage(normalize(ref_beam.imag))

        print("Ref x hologram: normalized / equalized")
        rhologram = ref_beam * hologram
        showimage(normalize(rhologram), equalize(rhologram))
        
        chirp = get_chirp(shape, distance)
        delilah = chirp * rhologram
        delilah_dft = get_shifted_dft(delilah)
        
        wfa = get_wrapped_formula_array(shape, distance)
        reconstructed = wfa * delilah_dft

        module = normalize(abs(reconstructed))
        print("abs(fresnel)")
        showimage(module)

        phase = angle2(reconstructed)
        print("phase(fresnel)")
        showimage(normalize(phase))
    return 0


if __name__ == "__main__":
    exit(main())
