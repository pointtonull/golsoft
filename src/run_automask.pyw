#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from automask import get_mask, get_circles, get_holed_window
from autopipe import showimage, red, blue
from scipy import misc
import Image as pil
import numpy as np
import sys

HOLE_LEN = 4


def kaiser14(length):
    return np.kaiser(length, 14)



def main():

    windows = {
        "Rectangular": np.ones,
        "Hanning": np.hanning,
        "Bartlett": np.hamming,
        "Blackman": np.bartlett,
        "Hamming": np.blackman,
        "Kaiser14": kaiser14,
    }

    if len(sys.argv) < 2:

        for w_name in windows:

            print("%s window:" % w_name)
            window = windows[w_name](500)
            window2d = get_mask((500, 500), window)
            showimage(window2d * 255)

    else:

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
