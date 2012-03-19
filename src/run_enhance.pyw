#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from autopipe import red, blue, showimage
from enhance import logscale, equalize, get_centered
from scipy import misc
import numpy as np
import sys


def main():
    images = [misc.imread(filename, True) for filename in sys.argv[1:]]
    if not images:
        lena = misc.imresize(misc.lena(), .5)
        images = [lena]

    for image in images:
        width, height = image.shape
        if max(width, height) > 600:
            print("Resizing...")
            prop = 600. / max(width, height)
            image = misc.imresize(image, prop)

        print("Original image:")
        showimage(image)

        print("With center of mass to center of image:")
        centered = get_centered(image)
        showimage(centered)

        print("With logscale:")
        showimage(logscale(image))

        print("With equalized histogram:")
        showimage(equalize(image))

if __name__ == "__main__":
    exit(main())
