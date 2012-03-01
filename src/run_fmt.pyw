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
from enhance import equalize, logscale, get_intensity
from fmt import get_logpolar, get_fmt, get_fmt_correlation
from itertools import product, combinations, permutations
from random import sample
from scipy import misc, ndimage
import sys



def main():
    images = [misc.imread(filename, True) for filename in sys.argv[1:]]
    if not images:
        lena = misc.imresize(misc.lena(), .5)
        images = [lena]
    angles = (0., 22.5, 45, 90., 135., 180., 270)
    scales = (.5, .75, 1.)
    translations = product((-15, 0, 15), (-15, 0, 15))

    transformations = []
    samples = sample(list(product(images, scales, angles, translations)), 5)
    for combination in samples:
        image, scale, angle, translation = combination
        image = misc.imrotate(image, angle)
        image = ndimage.shift(image, translation)
        image = misc.imresize(image, scale)
        transformations.append(image)

    for image1, image2 in combinations(transformations, 2):
        print("\nCompare images:")
        showimage(image1, image2)
        fmt1 = equalize(get_intensity(get_fmt(image1)))
        fmt2 = equalize(get_intensity(get_fmt(image2)))
        showimage(fmt1, fmt2)
        print(get_fmt_correlation(image1, image2))


if __name__ == "__main__":
    exit(main())
