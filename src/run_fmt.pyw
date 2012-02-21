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
import sys



def main():
    images = [misc.imread(filename, True) for filename in sys.argv[1:]]
    if not images:
        images = [misc.lena()]
    for image in images:
        rpi = misc.imrotate(image, 90)
        rtau = misc.imrotate(image, 180)
        for transform in (image, rpi, rtau):
            showimage(transform)
            fmt = get_fmt(transform)
            showimage(equalize(fmt.real ** 2 + fmt.imag ** 2))


if __name__ == "__main__":
    exit(main())
