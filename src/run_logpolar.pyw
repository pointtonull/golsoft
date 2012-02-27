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
from fmt import get_logpolar
from scipy import misc
from itertools import product
import sys



def main():
    images = [misc.imread(filename, True) for filename in sys.argv[1:]]
    if len(images) < 2:
        if not images:
            lena = misc.imresize(misc.lena(), .5)
            images = [lena]

    for image in images:
        showimage(image)
        logpolar = get_logpolar(image, 3)
        showimage(logpolar)

if __name__ == "__main__":
    exit(main())
