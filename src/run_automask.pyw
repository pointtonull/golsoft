#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from automask import get_mask, get_circles, get_holed_window
from autopipe import showimage, red, blue
from itertools import product
from scipy import misc
import Image as pil
import numpy as np
from numpy import kaiser
import sys
from enhance import equalize


def main():

    alone = [False]
    radius_factor = [1.5, 2, 2.5]
    softness = [2, 4, 8]
    hole_len = [4, 6]

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
            circles = sorted((get_circles(array)))
            array = equalize(array)
            print("    %s" % circles)
            for value, (crow, ccol), radius in circles:
                for alone, radius_factor, softness, hole_len in options:
                    scaled_radius = radius * radius_factor
                    print("        Center %d, %d, radius %d:" %
                        (ccol, crow, scaled_radius))
                    top, bottom = crow - radius, crow + scaled_radius
                    left, right = ccol - radius, ccol + scaled_radius
                    helprect = "X%d-%d Y%d-%d" % (left, right, top, bottom)
                    print("        %s" % helprect)
                    optsstring = "hole%d %sks%0.1f" % (
                        hole_len, "Alone " * alone, softness)
                    filename = "-".join((infilename, optsstring, helprect))
                    filename += ".tiff"
                    print(filename)
                    windowmaker = lambda x: kaiser(x, softness)
                    window = get_holed_window(windowmaker, scaled_radius,
                        hole_len)
                    mask = get_mask(array.shape, window, (crow, ccol))
                    image = pil.fromarray(np.float32(mask))
                    image.save(filename)
                    showimage(mask * 255)
                    showimage(mask * array)


if __name__ == "__main__":
    exit(main())
