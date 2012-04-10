#!/usr/bin/env python
#-*- coding: UTF-8 -*-

"""
Simple aplication of to show the Fourier transforms of images
"""

from autopipe import showimage
from dft import get_shifted_dft, get_shifted_idft
from dft import get_dft, get_idft
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
        showimage(image)

        dft_complex = get_shifted_dft(image)
        dft_intensity = equalize2(dft_complex)
        print("Intensity:")
        showimage(dft_intensity)

        rec_complex = get_shifted_idft(dft_complex)
        rec_image = normalize((rec_complex.real))
        print("Reconstructed:")
        showimage(rec_image)

if __name__ == "__main__":
    exit(main())
