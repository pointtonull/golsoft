#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from autopipe import showimage
from image import equalize
from fmt import get_fmt, get_fmt_correlation
from itertools import product, combinations
from random import sample
from scipy import misc, ndimage
import sys
import numpy as np

def main():
    images = [misc.imread(filename, True) for filename in sys.argv[1:]]
    if not images:
        lena = misc.imresize(misc.lena(), .5)
        images = [lena]
    angles = np.random.sample(10) * 360
    scales = (.25, .5, .75, 1.)
    translations = product((-15, 0, 15), (-15, 0, 15))

    transformations = []
    samples = sample(list(product(images, scales, angles, translations)), 20)
    for combination in samples:
        image, scale, angle, translation = combination
        image = misc.imrotate(image, angle)
        image = ndimage.shift(image, translation)
        image = misc.imresize(image, scale)
        transformations.append(image)

    print("Testing comparations")
    for image1 in transformations:
        print("\n------------------")
        print("Specimen / result")
        results = [(get_fmt_correlation(image1, image2, 0.3), image2)
            for image2 in images]
        total = sum((result[0][0] for result in results))
        best_score, best_match = max(results)
        showimage(image1, best_match)
        print("scores:")
        for result in sorted(results, reverse=True):
            (score, pos), image = result
            print("    %f, %s" % (score / total, str(pos)))


if __name__ == "__main__":
    exit(main())
