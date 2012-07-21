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
from pea import get_auto_mask, angle2
from ranges import frange
from unwrap import unwrap_wls, unwrap_qg


def main():
    """
    Le main rutine
    """
    images = [(filename, imread(filename, True))
        for filename in sys.argv[1:]]

    for filename, hologram in images:
        print(filename)
        showimage(equalize(hologram))

#        cos_alpha, cos_beta = calculate_director_cosines(hologram)
#        print("Cosines: %4.3f %4.3f" % (cos_alpha, cos_beta))

        distances = (-.05, )
        evenshapes = (False, )
        for distance, to_evenshape in zip(distances, evenshapes):
            print("Distance: %3.2f" % distance)
            if to_evenshape:
                image = evenshape(hologram, True)
            else:
                image = hologram

            shape = image.shape
            print("Shape: %d, %d" % shape)

#        print("Reference beam: normalized(imag == real)")
#        ref_beam = get_ref_beam(shape, cos_alpha, cos_beta)
#        showimage(normalize(ref_beam.imag))
#
#        print("Ref x hologram: normalized / equalized")
#        rhologram = ref_beam * hologram

            rhologram = image
            showimage(normalize(rhologram), equalize(rhologram))

            print("Spectrum:")
            spectrum = get_shifted_dft(rhologram)
            intensity = get_intensity(spectrum)
            showimage(equalize(intensity))

            print("Masked spectrum")
            softness = 0
            print("Mask softness; %3.2f" % softness)

            radious_scale = 1.5
            print("Radious scale; %3.2f" % radious_scale)

            zero_scale = 1.6
            mask, masked, centered = get_auto_mask(spectrum, softness,
                radious_scale, zero_scale)

            showimage(normalize(mask), normalize(equalize(masked)
                + equalize(centered)))

            print("To evenshape: %s" % to_evenshape)

            propagation_array = get_propagation_array(shape, distance)
            propagated = propagation_array * masked

            reconstructed = get_shifted_idft(propagated)
            module = normalize(abs(reconstructed))
            print("Module")
            showimage(module)
#            imwrite(module, "%s-module.jpg" % filename)
            phase = angle2(reconstructed)
            print("Phase")
            showimage(normalize(phase))
#            imwrite(phase, "%s-phase.jpg" % filename)
            unwrapped = unwrap_wls(phase)
            print("Least Squares Unwrapped Phase")
            showimage(normalize(unwrapped))

#            import statprof
#            statprof.start()
            
            unwrapped = unwrap_qg(phase, module)

#            statprof.stop()
#            statprof.display()

            print("Quality Guided Unwrapped Phase")
            showimage(normalize(unwrapped))
    return 0


if __name__ == "__main__":

    exit(main())
