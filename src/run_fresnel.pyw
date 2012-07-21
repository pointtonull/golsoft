# !/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Testing case
"""

from autopipe import showimage
from fresnel import angle2, get_chirp, get_wrapped_formula_array
import sys
from image import equalize, imread, normalize
from pea import calculate_director_cosines, get_ref_beam
from dft import get_shifted_dft
import numpy as np



def main():
    """
    Le main rutine
    """
    images = [(filename, imread(filename, True))
        for filename in sys.argv[1:]]
    
    for filename, hologram in images:
        print(filename)
        for porcion in (512, 256):
            hologram = hologram[:porcion, :porcion]
            shape = hologram.shape
            hologram = hologram - hologram.mean()
            showimage(equalize(hologram))
            
            distance = .16
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
