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
from itertools import product
from scipy import misc
import sys



def main():
    images = [misc.imread(filename, True) for filename in sys.argv[1:]]
    if len(images) < 2:
        if not images:
            lena = misc.imresize(misc.lena(), .5)
            images = [lena]

    for image in images:
        showimage(image)

        logpolar5 = get_logpolar(image, 5)
        recart5 = get_logpolar(logpolar5, 5, True)
        showimage(logpolar5, recart5)


if __name__ == "__main__":
    exit(main())
