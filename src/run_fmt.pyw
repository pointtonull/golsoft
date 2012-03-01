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
    #TODO: center image using center of mass of the image
    images = [misc.imread(filename, True) for filename in sys.argv[1:]]
    if len(images) < 2:
        if not images:
            lena = misc.imresize(misc.lena(), .5)
            images = [lena]
        angles = (0., 90., 180.)
        scales = (.5, 1.)

        for combination in product(images, scales, angles):
            image, scale, angle = combination
            print(scale, angle)
            image = misc.imrotate(image, angle)
            image = misc.imresize(image, scale)
            fmt = get_fmt(image)
            fmt_intensity = equalize(fmt.real ** 2 + fmt.imag ** 2)
            showimage(image, fmt_intensity)
    else:
        for image1, image2 in combinations(images, 2):
            print("\nCompare images:")
            showimage(image1, image2)
            get_correlation(image1, image2)


if __name__ == "__main__":
    exit(main())
