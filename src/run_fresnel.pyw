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

        distance = .11
        print("Distance: %3.2f" % distance)

#        print("Reference beam: normalized(imag == real)")
#        ref_beam = get_ref_beam(shape, cos_alpha, cos_beta)
#        showimage(normalize(ref_beam.imag))
#
#        print("Ref x hologram: normalized / equalized")
#        rhologram = ref_beam * hologram
        rhologram = hologram
        showimage(normalize(rhologram), equalize(rhologram))

        print("Spectrum:")
        spectrum = get_shifted_dft(rhologram)
        intensity = get_intensity(spectrum)
        showimage(equalize(intensity))

        print("Masked spectrum")
        softness = 1
        print("Mask softness; %3.2f" % softness)

        radious_scale = 1.8
        print("Radious scale; %3.2f" % radious_scale)
        
        tau = 6.283185307179586
        pi = 3.1415926589
        LAMBDA = 6.328e-07 # wave length
        DX = 8.39e-6
        DY = 8.46e-6
        K = tau / LAMBDA # wave number
        EPSILON = 1e-16


        def get_chirp(shape):
            maxrow = shape[0] / 2
            maxcol = shape[1] / 2
            minrow, mincol = -maxrow, -maxcol
            row, col = np.ogrid[minrow:maxrow:1., mincol:maxcol:1.]
            fresnel = exp(-1j * pi/(LAMBDA * distance) * (col ** 2 * DX ** 2
                + row ** 2 * DY ** 2))
            return fresnel

            
        def get_wrapped_formula_array(shape):
            maxrow = shape[0] / 2
            maxcol = shape[1] / 2
            minrow, mincol = -maxrow, -maxcol
            row, col = np.ogrid[minrow:maxrow:1., mincol:maxcol:1.]
            eta=LAMBDA*distance/(512*DX)
            nu=LAMBDA*distance/(512*DY)
            wfa = (1j / (LAMBDA * distance)) * exp(-1j * pi/LAMBDA
                *(2*distance+1/distance** (row ** 2 * eta ** 2 + col** 2 * nu ** 2)))
            return wfa
            

        cuttop = 0
        print("Cuttop: %2.2f" % cuttop)

        zero_scale = 1.6
        mask, masked = get_auto_mask(spectrum, softness, radious_scale,
            zero_scale, cuttop)

        print("Mask")
        showimage(normalize(mask), equalize(masked))

        masked_idft = get_shifted_idft(masked)

        chirp = get_chirp(shape)
        
        print("abs(chirp)")
        showimage(equalize(abs(chirp)))

#        delilah = chirp * masked_idft
        delilah = chirp * hologram
        delilah_dft = get_shifted_dft(delilah)

        print("intensity(delilah_dft)")
        showimage(equalize(delilah_dft))

        wfa = get_wrapped_formula_array(shape)
        print("intensity(wfa)")
        showimage(equalize(wfa))

        reconstructed = wfa * delilah_dft

        module = normalize(abs(reconstructed))
        print("abs(fresnel)")
        showimage(module)
#            imwrite(module, "%s-module.jpg" % filename)
        phase = angle2(reconstructed)
        print("phase(fresnel)")
        showimage(normalize(phase))
#            imwrite(phase, "%s-phase.jpg" % filename)
    return 0


if __name__ == "__main__":
    exit(main())
