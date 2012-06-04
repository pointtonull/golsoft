#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from autopipe import showimage
from image import imread, derotate, get_centered
from scipy.misc import imresize, imrotate, lena
from itertools import product
from random import sample
import sys


def main():
    images = [(filename, get_centered(imread(filename)))
        for filename in sys.argv[1:]]
    angles = range(360)

    if not images:
        image = imresize(lena(), .5)
        images = ["lena", image]

    samples = sample(list(product(images, angles)), 20)

    transformations = []
    for combination in samples:
        (filename, image), angle = combination
        image = imrotate(image, angle)
        transformations.append(image)

    print("Rotated / derotated image:")
    for image in transformations:
        derotated = derotate(image)
        showimage(image, derotated)

    return 0

if __name__ == "__main__":
    exit(main())
