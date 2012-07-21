#!/usr/bin/env python
#-*- coding: UTF-8 -*-

import sys

from scipy.misc import imresize, lena

from autopipe import showimage
from image import imread, normalize
from unwrap import unwrap_wls, unwrap_goldstein

def read(filename):
    try:
        return imread(filename)
    except:
        return None

def main():
    images = [(filename, read(filename))
        for filename in sys.argv[1:]]

    unwrappers = (
        ("WLS", unwrap_wls),
        ("Goldstein", unwrap_goldstein),
    )

    if not images:
        image = imresize(lena(), .5)
        images = ["lena", image]

    for filename, image in images:
        if image is not None:
            for name, unwrapper in unwrappers:
                print("Unwrapping %s with %s:" % (filename, name))
                unwrapped = unwrapper(image)
                showimage(image, normalize(unwrapped))

    return 0

if __name__ == "__main__":
    exit(main())
