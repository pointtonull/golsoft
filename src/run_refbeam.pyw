#!/usr/bin/env python
#-*- coding: UTF-8 -*-

"""
This is a simple implementation of the Fourier-Mellin Trasform
image = [m, n]
ft_magnitude = |fft(image)|
lp_ft_magnitude = logpolar(ft_magnitude)
fmt = fft(lp_ft_magnitude)
"""

import sys

from numpy import hstack, sin, cos
from scipy import misc

from autopipe import showimage
from dft import get_shifted_dft
from image import imread, normalize
from pea import get_refbeam, calculate_director_cosines, PEA


def main():
    images = [(filename, imread(filename, True)) for filename in sys.argv[1:]]
    if len(images) < 2:
        if not images:
            lena = misc.imresize(misc.lena(), .5)
            images = [("lena", lena)]

    for filename, image in images:
        print("Original image: %s" % filename)
        image = normalize(image)
        methods = (
            ("Direct method", calculate_director_cosines),
        )

        for description, function in methods:

            cos_alpha, cos_beta = function(get_shifted_dft(image), 1, (1, 1))
            print("\n%s\n" % description)
            print(cos_alpha, cos_beta)

            pea = PEA(filename)
            phase = pea.phase
            ref_beam = get_refbeam(image.shape, cos_alpha, cos_beta, 1, (1, 1))
            showimage(hstack((normalize(image), normalize(ref_beam.imag))))
            showimage(hstack(((ref_beam.imag +  cos(phase)) % 1, ref_beam.imag)))

    return 0


if __name__ == "__main__":
    exit(main())
