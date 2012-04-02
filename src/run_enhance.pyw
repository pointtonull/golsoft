#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from autopipe import red, blue, showimage
from image import logscale, equalize, get_centered, imread
import numpy as np
import sys
from scipy.misc import imresize


def main():
    images = [(filename, imread(filename)) for filename in sys.argv[1:]]
    if not images:
        lena = imresize(misc.lena(), .5)
        images = [lena]

    for filename, image in images:
        print(filename)
        width, height = image.shape
        if max(width, height) > 600:
            print("Resizing...")
            prop = 600. / max(width, height)
            image = imresize(image, prop)

#        print("Original image:")
#        showimage(image)

#        print("With center of mass to center of image:")
#        centered = get_centered(image)
#        showimage(centered)

        print("With logscale:")
        showimage(logscale(image))

#        print("With equalized histogram:")
#        showimage(equalize(image))
    return 0

if __name__ == "__main__":
    exit(main())
