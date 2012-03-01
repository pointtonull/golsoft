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
from fmt import get_logpolar, cv_logpolar
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

#        logpolar0 = get_logpolar(image, 0)
#        recart0 = get_logpolar(logpolar0, 0, True)
#        showimage(logpolar0, recart0)

#        logpolar1 = get_logpolar(image, 1)
#        recart1 = get_logpolar(image, 1, True)
#        showimage(logpolar1, recart1)

#        logpolar2 = get_logpolar(image, 2)
#        recart2 = get_logpolar(image, 2, True)
#        showimage(logpolar2, recart2)

#        logpolar3 = get_logpolar(image, 3)
#        recart3 = get_logpolar(image, 3, True)
#        showimage(logpolar3, recart3)

#        logpolar4 = get_logpolar(image, 4)
#        recart4 = get_logpolar(image, 4, True)
#        showimage(logpolar4, recart4)

        logpolar5 = get_logpolar(image, 5)
        recart5 = get_logpolar(logpolar5, 5, True)
        showimage(logpolar5, recart5)


if __name__ == "__main__":
    exit(main())
