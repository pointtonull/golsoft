#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from automask import get_mask, get_circles, get_holed_window
from autopipe import showimage, red, blue
from itertools import product
from scipy import misc
import Image as pil
import numpy as np
import sys


def main():

    alone = [True, False]
    radius_factor = [1., 1.25, 1.5, 1.75, 2]
    softness = [0, 2, 5, 6, 8.6]
    hole_len = [0, 4]

    if len(sys.argv) < 2:
        options = product(solitudiness, softness, hole_len)
        for solitudiness, softness, hole_len in options:
            print("%so %s %0.2fs") % ("alone" * alone, softness)
            window = np.kaiser(250, softness)
            window2d = get_mask((250, 250), window)
            showimage(window2d * 255)
    else:
        options = product(alone, radius_factor, softness, hole_len)
        for infilename in sys.argv[1:]:
            print("File: %s" % infilename)
            array = misc.imread(infilename)
            if "." in infilename:
                infilename = "".join(infilename.split(".")[:-1])
            circles = sorted((get_circles(array)))
            for value, (crow, ccol), radius in circles:
                print("  Center %d, %d, radius %d:" % (ccol, crow, radius))
                top, bottom = crow - radius, crow + radius
                left, right = ccol - radius, ccol + radius
                helprect = "X%d-%d Y%d-%d" % (left, right, top, bottom)
                print("    %s" % helprect)
                for w_name, w_function in windows.iteritems():
                    plusholes = {
                        "o" + w_name: lambda length:get_holed_window(w_function,
                            length, HOLE_LEN),
                        w_name: w_function,
                    }
                    for w_name, w_func in plusholes.iteritems():
                        filename = "-".join((infilename, w_name, helprect))
                        filename += ".tiff"
                        print("      %s" % w_name)
                        window = w_func(radius * 2)
                        mask = get_mask((512, 512), window, (crow, ccol))
                        image = pil.fromarray(np.float32(mask))
                        image.save(filename)
                        showimage(mask * 255)


if __name__ == "__main__":
    exit(main())
