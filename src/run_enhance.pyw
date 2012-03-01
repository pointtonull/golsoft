#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from autopipe import red, blue, showimage
from enhance import logscale, equalize
import Image
import sys


def main():
    images = [misc.imread(filename, True) for filename in sys.argv[1:]]
    if not images:
        lena = misc.imresize(misc.lena(), .5)
        images = [lena]

    for image in images:
        width, height = image.size
        if max(width, height) > 600:
            print("Resizing...")
            prop = max(width, height) / 600.
            image = image.resize((int(width / prop), int(height / prop)), 1)

        print("Original image:")
        showimage(image)

        print("With logscale:")
        showimage(logscale(image))

        print("With equalized histogram:")
        showimage(Image.fromarray(equalize(image)))

if __name__ == "__main__":
    exit(main())
