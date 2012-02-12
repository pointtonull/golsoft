#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from autopipe import showimage, red, blue
import numpy as np
from scipy.ndimage import maximum_filter1d, minimum_filter1d
import Image as pil
import sys
import pylab

def get_maxs(array):
    xsums = np.array([array[:, x].sum()
        for x in range(array.shape[-1])])
    ddiffs = np.diff(xsums)
    localmaxs = maximum_filter1d(ddiffs, 10)
    localmins = minimum_filter1d(localmaxs, 10)
    pylab.plot(xsums)
    pylab.plot(localmins)
    print(localmins.argsort()[-3:])
    pylab.show()

def main():
    for filename in sys.argv[1:]:
        image = pil.open(filename)
        array = np.asarray(image)
        get_maxs(array)

if __name__ == "__main__":
    exit(main())
