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
from fmt import get_logpolar, get_fmt
from scipy import misc
from itertools import product, combinations
from scipy.ndimage import shift
import sys



def main():
    images = [misc.imread(filename, True) for filename in sys.argv[1:]]
    if len(images) < 2:
        if not images:
            images = [misc.lena()]
        angles = (45, 90, 180)
        scales = (.5, 1., 1.5)
        centers = list(combinations((-15, 0, 15), 2))
        for combination in product(images, angles, scales, centers):
            image, angle, scale, center = combination
            image = misc.imrotate(image, angle)
            image = misc.imresize(image, scale)
            image = shift(image, center)
            showimage(image)
            fmt = get_fmt(image)
            showimage(equalize(fmt.real ** 2 + fmt.imag ** 2))
    else:
        print("Compare images")


if __name__ == "__main__":
    exit(main())
