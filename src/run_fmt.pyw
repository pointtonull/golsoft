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
import sys



def main():
    images = [misc.imread(filename, True) for filename in sys.argv[1:]]
    if not images:
        images = [misc.lena()]
    for image in images:
        logpolar = get_logpolar(image)
#        showimage(image)
        showimage(logpolar)


if __name__ == "__main__":
    exit(main())
