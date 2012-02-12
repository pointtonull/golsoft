#!/usr/bin/env python
#-*- coding: UTF-8 -*-

"""
This is a simple implementation of the Fourier-Mellin Trasform
image = [m, n]
ft_magnitude = |fft(image)|
lp_ft_magnitude = logpolar(ft_magnitude)
fmt = fft(lp_ft_magnitude)
"""

from autopipe import showimage
from fmt import logpolar, asarray
import Image
import sys


def main():
    filename = sys.argv[1]
    print("Openning image file: %s" % filename)
    image = Image.open(filename)
    showimage(image)
    print("Transforming to log polar")
    image = Image.fromarray(logpolar(asarray(image), 1))
    print("Displaying")
    showimage(image)

if __name__ == "__main__":
    exit(main())
