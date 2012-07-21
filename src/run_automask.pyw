#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from automask import get_mask, get_circles, get_holed_window
from autopipe import showimage, red, blue
from dft import get_shifted_dft
from image import equalize, logscale, imread
from itertools import product
from numpy import kaiser
from scipy import misc
import numpy as np
import sys


def main():

    alone = [True, False]
    radius_factor = [2,  3]
    softness = [2, 4]
    hole_len = [4, 8]

    images = [(filename, imread(filename, True)) for filename in sys.argv[1:]]
    for filename, image in images:
        print("File: %s" % filename)
        image = equalize(get_shifted_dft(image))
        circles = get_circles(image, 5)
        l_pnoise, l_order, c_order, r_order, r_pnoise = circles

        windowmaker = lambda x: kaiser(int(round(x)), softness)
        window = get_holed_window(windowmaker, l_order[2] * r_scale,
            hole_len)
        mask = get_mask(array.shape, window, l_order[1])
        if alone:
            window = windowmaker(l_pnoise[2] * r_scale)
            mask *= 1 - get_mask(array.shape, window, l_pnoise[1])
            window = windowmaker(c_order[2] * r_scale)
            mask *= 1 - get_mask(array.shape, window, c_order[1])
        showimage(mask * 255)
        image = pil.fromarray(np.float32(mask))
        image.save(filename + " Left.tiff")

        window = get_holed_window(windowmaker, r_order[2] * r_scale,
            hole_len)
        mask = get_mask(array.shape, window, r_order[1])
        if alone:
            window = windowmaker(r_pnoise[2] * r_scale)
            mask *= 1 - get_mask(array.shape, window, r_pnoise[1])
            window = windowmaker(c_order[2] * r_scale)
            mask *= 1 - get_mask(array.shape, window, c_order[1])
        showimage(mask * 255)
        image = pil.fromarray(np.float32(mask))
    return 0

if __name__ == "__main__":
    exit(main())
