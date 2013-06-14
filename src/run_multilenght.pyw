#!/usr/bin/env python
#-*- coding: UTF-8 -*-

import sys

from scipy import misc
import numpy as np
from numpy import hstack

from autopipe import showimage
from image import normalize, equalize, imread, imwrite
from automask import get_auto_mask
import dft
import pea

tau = 2 * np. pi


def main():
    images = [(filename, imread(filename, True)) for filename in sys.argv[1:]]
    if len(images) < 2:
        if not images:
            lena = misc.imresize(misc.lena(), .5)
            images = [("lena", lena)]

    for filename, depth in images:
        print("Depth map: %s" % filename)
        print("shape: %d, %d" % depth.shape)

        inclined_plane = pea.get_inclined_plane(depth.shape, 0.90, 0.90, 20, (1, 1))
        depth = inclined_plane + depth / 2. ** 6
        
        holo1 = pea.get_fringes(depth, 2)
        holo2 = pea.get_fringes(depth, 1.5)
        showimage(hstack((normalize(depth - inclined_plane * 0.9), normalize(holo1),
            normalize(holo2))))

        spectrum1 = dft.get_shifted_dft(holo1)
        mask, masked, centered = get_auto_mask(spectrum1, centered=False)
        showimage(equalize(masked) + equalize(spectrum1), equalize(spectrum1))
        r1 = p.reconstructed
        p.image = holo2
        r2 = p.reconstructed
        rs = r1 + r2
        holo3 = pea.get_module(rs)
        p.image = holo3
        showimage(p.phase_corrected)

    return 0


if __name__ == "__main__":
    exit(main())
