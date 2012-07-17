
# !/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Testing case
"""

import sys

from numpy import abs, arctan2
import numpy as np

from autopipe import showimage
from dft import get_shifted_dft, get_shifted_idft
from image import equalize, imread, normalize, imwrite, get_intensity
from image import evenshape
from pea import calculate_director_cosines, get_ref_beam, get_propagation_array
from pea import get_auto_mask
from ranges import frange
from unwrap import unwrap_wls, unwrap_qg

def angle2(array):
    return arctan2(array.real, array.imag)

def main():
    """
    Le main rutine
    """
    images = [(filename, imread(filename, True))
        for filename in sys.argv[1:]]

    for filename, hologram in images:
        print(filename)
        showimage(hologram)
        shape = hologram.shape
        spectrum = get_shifted_dft(hologram)

        distances = (-11, -.05, 0, .05, .11)
        for distance in distances:
            print("Distance: %3.2f" % distance)

            propagation_array = get_propagation_array(shape, distance)
            propagated = propagation_array * spectrum

            reconstructed = get_shifted_idft(propagated)
            module = normalize(abs(reconstructed))
            showimage(module, normalize(angle2(propagated)))
    return 0


if __name__ == "__main__":

    exit(main())
