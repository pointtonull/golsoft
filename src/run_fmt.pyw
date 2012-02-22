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
from enhance import equalize
from fmt import get_logpolar, get_fmt, get_correlation
from scipy import misc
from itertools import product, combinations
from scipy.ndimage import shift
import sys



def main():
    images = [misc.imread(filename, True) for filename in sys.argv[1:]]
    if len(images) < 2:
        if not images:
            lena = misc.imresize(misc.lena(), .5)
            images = [lena]
        angles = (0., 90., 180.)
        scales = (.5, 1., 1.5)
        scales = (.5, 1.)
        centers = list(combinations((-25, 0, 25), 2))

        for combination in product(images, centers, scales, angles):
            image, center, scale, angle = combination
            print(scale, center, angle)
            image = shift(image, center)
            image = misc.imrotate(image, angle)
            image = misc.imresize(image, scale)
            fmt = get_fmt(image)
            fmt_intensity = equalize(fmt.real ** 2 + fmt.imag ** 2)
            showimage(image, fmt_intensity)
    else:
        for image1, image2 in combinations(images, 2):
            print("\nCompare images:")
            showimage(image1, image2)
            print("Correlation: %0.5f" % get_correlation(image1, image2))


if __name__ == "__main__":
    exit(main())
