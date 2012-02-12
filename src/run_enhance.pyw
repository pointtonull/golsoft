#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from autopipe import red, blue, showimage
from enhance import logscale, equalize, autocontrast
import Image
import sys


def main():
    image = Image.open(sys.argv[1])
    width, height = image.size
    if max(width, height) > 600:
        print("Resizing...")
        prop = max(width, height) / 600.
        image = image.resize((int(width / prop), int(height / prop)), 1)

    print("Original image:")
    showimage(image)

    print("With logscale:")
    showimage(logscale(image))

    print("With auto-contrast:")
    showimage(autocontrast(image))

    print("With equalized histogram:")
    showimage(equalize(image))

if __name__ == "__main__":
    exit(main())
