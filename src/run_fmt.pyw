#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from autopipe import showimage
from enhance import equalize, logscale
from fmt import get_logpolar, get_fmt, get_fmt_correlation
from itertools import product, combinations, permutations
from random import sample
from scipy import misc, ndimage
import sys


def main():
    images = [misc.imread(filename, True) for filename in sys.argv[1:]]
    if not images:
        lena = misc.imresize(misc.lena(), .5)
        images = [lena]
    angles = (-15., -10., -5., 0., 5., 10., 15.)
    scales = (.75, 1., 1.25)
    translations = product((-15, 0, 15), (-15, 0, 15))

    transformations = []
    samples = sample(list(product(images, scales, angles, translations)), 3)
    for combination in samples:
        image, scale, angle, translation = combination
        image = misc.imrotate(image, angle)
        image = ndimage.shift(image, translation)
        image = misc.imresize(image, scale)
        transformations.append(image)

    for image1, image2 in combinations(transformations, 2):
        print("\nCompare images:")
        showimage(image1, image2)
        fmt1 = equalize(get_fmt(image1))
        fmt2 = equalize(get_fmt(image2))
        showimage(fmt1, fmt2)
        print(get_fmt_correlation(image1, image2))


if __name__ == "__main__":
    exit(main())
