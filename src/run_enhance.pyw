#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from autopipe import showimage
from image import equalize, imread, get_centered, derotate
from scipy.misc import imresize
from itertools import product
from random import sample
import sys


def main():
    images = [(filename, imread(filename)) for filename in sys.argv[1:]]
    angles = (-15., -10., -5., 0., 5., 10., 15.)

    if not images:
        lena = imresize(misc.lena(), .5)
        images = [lena]

    samples = sample(list(product(images, angles)), 3)

    transformations = []
    for combination in samples:
        image, angle = combination
        image = misc.imrotate(image, angle)
        transformations.append(image)

    for image in transformations:
        print("\nCompare images:")
        showimage(image)
        derotated = derotate(image)
        showimage(derotated)

    return 0

if __name__ == "__main__":
    exit(main())
