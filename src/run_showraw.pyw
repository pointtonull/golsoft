#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from autopipe import showimage
from image import imread
from scipy.misc import imresize, lena
import sys


def main():
    images = [(filename, imread(filename))
        for filename in sys.argv[1:]]

    if not images:
        image = imresize(lena(), .5)
        images = ["lena", image]

    for filename, image in images:
        print(filename)
        showimage(image)

    return 0

if __name__ == "__main__":
    exit(main())
