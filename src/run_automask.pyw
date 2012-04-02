#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from automask import get_mask, get_circles, get_holed_window
from autopipe import showimage, red, blue
from image import equalize, logscale
from itertools import product
from numpy import kaiser
from scipy import misc
import Image as pil
import numpy as np
import sys


def main():

    alone = [True, False]
    radius_factor = [2,  3]
    softness = [2, 4]
    hole_len = [4, 8]

    if len(sys.argv) < 2:
        print("Invoqued without arguments, we will see some windows!\n")
        options = product(alone, softness, hole_len)
        length = 250
        for alone, softness, hole_len in options:
            print("hole%d %sks%0.1f") % (
                hole_len,
                "Alone " * alone,
                softness)
            windowmaker = lambda x: kaiser(x, softness)
            window = get_holed_window(windowmaker, length, hole_len)
            window2d = get_mask((length, length), window)
            showimage(window2d * 255)
    else:
        inputfiles = sys.argv[1:]
        print("Invoqued with %d arguments, trying to mask them..\n" %
            len(inputfiles))
        options = tuple(product(alone, radius_factor, softness, hole_len))
        for infilename in inputfiles:
            print("File: %s" % infilename)
            array = misc.imread(infilename)
            if "." in infilename:
                infilename = "".join(infilename.split(".")[:-1])
            circles = get_circles(array, 5)
            l_pnoise, l_order, c_order, r_order, r_pnoise = circles
            array = equalize(array)
            for alone, r_scale, softness, hole_len in options:
                optsstring = "Hole%d %sKS%0.1f" % (
                    hole_len, "Alone " * alone, softness)
                filename = "-".join((infilename, optsstring))
                print(filename)
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
                image.save(filename + " Right.tiff")
    return 0

if __name__ == "__main__":
    exit(main())
