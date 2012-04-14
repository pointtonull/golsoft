# !/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Testing case
"""

from autopipe import showimage
from image import equalize, imread, normalize
from pea import get_pea, calculate_director_cosines, guess_director_cosines
from numpy import abs, arctan2
from ranges import frange
import sys


def main():
    """
    Le main rutine
    """
    images = [(filename, imread(filename, True)) for filename in sys.argv[1:]]
    for filename, image in images:
        print(filename)
        showimage(equalize(image))
        cos_alpha, cos_beta = calculate_director_cosines(image)
        cos_alpha, cos_beta = guess_director_cosines(image)
#        cos_alpha, cos_beta = .98, .99
        for distance in frange(.0375, .005125, 1):
            print(distance)
            pea = get_pea(image, distance, cos_alpha, cos_beta)
            showimage(equalize(abs(pea)))
            showimage(normalize(arctan2(pea.real, pea.imag)))
    return 0


if __name__ == "__main__":
    exit(main())