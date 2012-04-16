# !/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Testing case
"""

from autopipe import showimage
from image import equalize, imread, normalize, imwrite
from pea import calculate_director_cosines, get_ref_beam, get_propagation_array
from pea import apply_mask
from numpy import abs, arctan2
from dft import get_shifted_dft, get_shifted_idft
import numpy as np
from ranges import frange
import sys

def angle2(array):
    return arctan2(array.real, array.imag)

def main():
    """
    Le main rutine
    """
    images = [(filename, imread(filename, True)) for filename in sys.argv[1:]]
    for filename, hologram in images:
        print(filename)
        shape = hologram.shape
        showimage(equalize(hologram))

        cos_alpha, cos_beta = calculate_director_cosines(hologram)
        print("Cosines: %4.3f %4.3f" % (cos_alpha, cos_beta))

        distance = .0375
        print("Distance: %3.2f" % distance)

        softness = 0
        print("Mask softness; %3.2f" % softness)

        radious_scale = 1.
        print("Radious scale; %3.2f" % radious_scale)

        print("Reference beam: normalized(imag == real)")
        ref_beam = get_ref_beam(shape, cos_alpha, cos_beta)
        showimage(normalize(ref_beam.imag))

        print("Ref x hologram: normalized / equalized")
        rhologram = ref_beam * hologram
        showimage(normalize(rhologram), equalize(rhologram))

        print("Spectrum:")
        spectrum = get_shifted_dft(rhologram)
        showimage(equalize(spectrum))

        print("Masked spectrum")
        for cuttop in frange(.5, .025, 1):
            print("Cuttop: %2.2f" % cuttop)
            masked_spectrum = apply_mask(spectrum, softness=softness,
                radious_scale=radious_scale, cuttop=cuttop)
            showimage(equalize(masked_spectrum))

            propagation_array = get_propagation_array(shape, distance)
            propagated = propagation_array * masked_spectrum

            reconstructed = get_shifted_idft(propagated)
            module = normalize(abs(reconstructed))
            showimage(module)
            imwrite(module, "%s-module.jpg" % (filename, cuttop))
            phase = angle2(reconstructed)
            showimage(normalize(phase))
            imwrite(phase, "%s-phase.jpg" % (filename, cuttop))
    return 0


if __name__ == "__main__":
    exit(main())
