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

    if not images:
        image = imresize(lena(), .5)
        images = ["lena", image]

    for image in images:
        showimage(image)

    return 0

if __name__ == "__main__":
    exit(main())
