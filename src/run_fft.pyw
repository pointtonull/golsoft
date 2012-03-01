#!/usr/bin/env python
#-*- coding: UTF-8 -*-

"""
This is a simple implementation of the Fourier-Mellin Trasform
image = [m, n]
ft_magnitude = |fft(image)|
lp_ft_magnitude = logpolar(ft_magnitude)
fmt = fft(lp_ft_magnitude)
"""

from autopipe import showimage, blue, red, green
from enhance import equalize, get_intensity
from fmt import get_shiftedfft
from scipy import misc, ndimage
import numpy as np
import sys



def main():
    images = [misc.imread(filename, True) for filename in sys.argv[1:]]
    if not images:
        lena = misc.imresize(misc.lena(), .5)
        images = [lena]

    for image in images:
        fft_complex = get_shiftedfft(image)
        fft_intensity = get_intensity(fft_complex)
        fft_image = equalize(fft_intensity)
        showimage(image, fft_image)

if __name__ == "__main__":
    exit(main())
