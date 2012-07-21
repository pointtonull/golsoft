#!/usr/bin/env python
#-*- coding: UTF-8 -*-

"""
This is a simple implementation of the Fourier-Mellin Trasform
image = [m, n]
ft_magnitude = |fft(image)|
lp_ft_magnitude = logpolar(ft_magnitude)
fmt = fft(lp_ft_magnitude)
"""

from autopipe import showimage
from image import imread, equalize, normalize, logscale
from scipy import misc
from pea import guess_director_cosines, get_ref_beam, calculate_director_cosines
from dft import get_shifted_dft
import sys



def main():
    images = [(filename, imread(filename, True)) for filename in sys.argv[1:]]
    if len(images) < 2:
        if not images:
            lena = misc.imresize(misc.lena(), .5)
            images = [("lena", lena)]

    for filename, image in images:
        print("Original image: %s" % filename)
        image = normalize(image)
        image_dft = get_shifted_dft(image)
        image_dft = equalize(image_dft)
        methods = (
            ("Naïve aproximation", guess_director_cosines),
            ("Direct method", calculate_director_cosines),
        )

        for description, function in methods:

            cos_alpha, cos_beta = function(image)
            print("\n%s\n" % description)
            print(cos_alpha, cos_beta)

            ref_beam = get_ref_beam(image.shape, cos_alpha, cos_beta)
            showimage(equalize(image), equalize(ref_beam.real))
            showimage(equalize(ref_beam.real), equalize(image))

            ref_beam_dft = get_shifted_dft(ref_beam.real)
            ref_beam_dft_cmp = 1 * equalize(ref_beam_dft)
            ref_beam_dft_cmp += 0 * logscale(ref_beam_dft)
            ref_beam_dft = normalize(ref_beam_dft_cmp)

            showimage(image_dft, ref_beam_dft)
            showimage(ref_beam_dft, image_dft)

    return 0


if __name__ == "__main__":
    exit(main())
