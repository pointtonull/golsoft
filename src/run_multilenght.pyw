#!/usr/bin/env python
#-*- coding: UTF-8 -*-

import sys

from scipy import misc
import numpy as np
from numpy import hstack

from automask import get_auto_mask
from autopipe import showimage
from image import normalize, equalize, imread, imwrite, get_centered, limit_size
from image import subtract
import dft
import pea
import unwrap

tau = 2 * np.pi


def main():
    images = [(filename, limit_size(imread(filename), 0.25))
        for filename in sys.argv[1:]]
    if not images:
        lena = misc.imresize(misc.lena(), .5)
        images = [("lena", lena)]

    if len(images) == 1:
        filename, depth = images[0]
        print("Depth map: %s" % filename)
        print("shape: %d, %d" % depth.shape)

        inclined_plane = pea.get_inclined_plane(depth.shape, 0.90, 0.90, 20,
            (1, 1))
        depth = inclined_plane + depth / 2. ** 6
        
        holo1 = normalize(pea.get_fringes(depth, 2.5))
        holo2 = normalize(pea.get_fringes(depth, 1.5))


    if len(images) == 2:
        
        file1, holo1 = images.pop()
        file2, holo2 = images.pop()

    p = pea.PEA()
    rows, cols = holo1.shape
    new_center = (rows * 3/4., cols * 3/4.)

    spectrum1 = dft.get_shifted_dft(holo1)
    mask, masked1, centered = get_auto_mask(spectrum1, centered=True)
    masked1 = get_centered(masked1, new_center, reverse=True)
    rec1 = dft.get_shifted_idft(masked1)
    rec1 /= pea.get_module(rec1)

    p.image = holo1
    phase1 = p.phase_corrected
    module1 = p.module

    spectrum2 = dft.get_shifted_dft(holo2)
    mask, masked2, centered = get_auto_mask(spectrum2, centered=True)
    masked1 = get_centered(masked2, new_center, reverse=True)
    rec2 = dft.get_shifted_idft(masked2)
    rec2 /= pea.get_module(rec2)

    p.image = holo2
    phase2 = p.phase_corrected
    module2 = p.module

    sintetic_holo = normalize(pea.get_module(rec1 + rec2))

    p.image_holo = sintetic_holo
#        showimage(hstack((normalize(p.image_holo), normalize(p.image_obj), normalize(p.image))))
    p.mask_zero_scale = 1.25
    p.mask_order_scale = 1
    try:
        spectrums = hstack((
            equalize(spectrum1),
            equalize(spectrum2),
            equalize(p.centered_spectrum) * 0.6 + equalize(p.spectrum_masked) * 0.4
        ))
    except:
        showimage(equalize(p.spectrum))
        raise

    holos = hstack((normalize(holo1), normalize(holo2), normalize(p.image)))
    phase3 = p.phase_corrected
    phases = hstack((normalize(phase1), normalize(phase2), normalize(phase3)))
    modules = normalize(hstack((module1, module2, p.module)))
    showimage(np.vstack((holos, spectrums, phases, modules)))

    congruent = unwrap.make_congruent(phase1, phase3)
    showimage(hstack((normalize(phase3), normalize(-congruent))))

#        p.unwrapper = unwrap.unwrap_qg
#        showimage(hstack((normalize(phase3), normalize(p.unwrapped_phase))))

    coarse = subtract(phase1, phase2)
    coarse = phase1 - phase2
    coarse[coarse < 0] += tau
    ignore, coarse = pea.align_phase(coarse)

    showimage(hstack((phase3, coarse)))
    return 0


if __name__ == "__main__":
    exit(main())
