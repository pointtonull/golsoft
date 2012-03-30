#!/usr/bin/env python
#-*- coding: UTF-8 -*-

"""
This is a simple implementation of the Fourier-Mellin Trasform
image = [m, n]
ft_magnitude = |fft(image)|
lp_ft_magnitude = logpolar(ft_magnitude)
fmt = fft(lp_ft_magnitude)
"""

from automask import get_circles, get_holed_window, get_mask
from autopipe import showimage, blue, red, green
from fmt import get_shiftedfft
from image import equalize, get_intensity, logscale, get_centered
from scipy import misc, ndimage
import numpy as np
import sys


MASK_SOFTNESS = 2
MASK_R_SCALE = 3

def apply_mask(array):
    shape = array.shape
    showimage(equalize(array), equalize(np.angle(array)), logscale(array.real))
    windowmaker = lambda x: np.kaiser(x, MASK_SOFTNESS)
    circles = sorted(get_circles(get_intensity(array), 3))
    virtual_order, real_order, zero_order = circles

#    for circle in circles:
#        print circle
#        centered = get_centered(array, circle[1])
#        showimage(equalize(centered))

    centered = get_centered(array, real_order[1])
    showimage(equalize(centered))
    print real_order[2]
    window = get_holed_window(windowmaker, real_order[2] * MASK_R_SCALE, 10)
    print shape, window.shape
    mask = get_mask(shape, window)

#    window = windowmaker(zero_order[2] * MASK_R_SCALE)
#    mask *= 1 - get_mask(shape, window)

#    showimage(logscale(mask * 255))
    masked = mask * centered
#    showimage(equalize(masked))
#    showimage(logscale(masked))
    return masked

def main():
    images = [misc.imread(filename, True) for filename in sys.argv[1:]]
    if not images:
        lena = misc.imresize(misc.lena(), .5)
        images = [lena]

    for image in images:
        fft_complex = get_shiftedfft(image)
        fft_intensity = equalize(fft_complex)
        showimage(image, fft_intensity)
        apply_mask(fft_complex)

if __name__ == "__main__":
    exit(main())
