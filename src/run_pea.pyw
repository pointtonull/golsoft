# !/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Testing case
"""

from autopipe import showimage
from image import equalize, imread, normalize
from pea import get_pea, guess_directors_angles
from numpy import angle
import sys


def frange(mean, radius, amount):
    """
    Aid for binary search of floats values
    """
    if amount == 1:
        start = mean
        step = 0.
    else:
        start = mean - radius
        step = (radius * 2) / (amount - 1.)

    for stepn in xrange(amount):
        yield start + stepn * step


def main():
    """
    Le main rutine
    """
    images = [(filename, imread(filename, True)) for filename in sys.argv[1:]]
    for filename, image in images:
        print(filename)
        showimage(equalize(image))
        for distance in frange(.0, .05, 2):
            print(distance)
            alpha, beta = guess_director_angles(image)
            pea = get_pea(image, distance, alpha, beta)
            showimage(equalize(pea))
            showimage(normalize(angle(pea)))
    return 0


if __name__ == "__main__":
    exit(main())
