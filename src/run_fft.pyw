#!/usr/bin/env python
#-*- coding: UTF-8 -*-

"""
Simple aplication of to show the Fourier transforms of images
"""

from autopipe import showimage
from fmt import get_shiftedfft
from image import equalize, imread
from scipy import misc
import sys


def main():
    """
    Le main rutine
    """
    images = [(filename, imread(filename)) for filename in sys.argv[1:]]
    if not images:
        lena = misc.imresize(misc.lena(), .5)
        images = [("lena", lena)]

    for filename, image in images:
        print("Original %s:" % filename)
        fft_complex = get_shiftedfft(image)
        fft_intensity = equalize(fft_complex)
        showimage(image, fft_intensity)

if __name__ == "__main__":
    exit(main())
