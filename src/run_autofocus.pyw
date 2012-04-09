#!/usr/bin/env python
#-*- coding: UTF-8 -*-

"""
This is a simple implementation of the Fourier-Mellin Trasform
image = [m, n]
ft_magnitude = |fft(image)|
lp_ft_magnitude = logpolar(ft_magnitude)
fmt = fft(lp_ft_magnitude)
"""

from StringIO import StringIO
from autopipe import showimage
from fmt import get_shiftedfft, get_ifft
from image import imread, equalize, normalize, logscale, get_intensity
from pea import apply_mask
from pea import guess_director_angles, get_ref_beam, get_propagation_array
from scipy import misc
import Image
from ranges import frange
import matplotlib.pyplot as plt
import numpy as np
import sys


def fig2image(figure):
    fileo = StringIO()
    figure.savefig(fileo)
    fileo.seek(0)
    image = Image.open(fileo)
    return image


def main():
    images = [(filename, imread(filename, True)) for filename in sys.argv[1:]]
    if len(images) < 2:
        if not images:
            lena = misc.imresize(misc.lena(), .5)
            images = [("lena", lena)]

    figure = plt.figure()
    for filename, image in images:

        print("Original image: %s" % filename)
        shape = image.shape
        showimage(image)

        alpha, beta = guess_director_angles(image)
        ref_beam = get_ref_beam(shape, alpha, beta)
        rhologram = ref_beam * image

        spectrum = get_shiftedfft(rhologram)
        masked = apply_mask(spectrum)
#        print(masked[0:4, 0:4])
#        showimage(normalize(masked[0:4, 0:4]))


        for distance in frange(-0.0359375
            , 2**-1, 11):
            print(distance)
            propagation_array = get_propagation_array(shape, distance)
#            print(propagation_array[0:4, 0:4])
#            showimage(normalize(masked[0:4, 0:4]))
#            showimage(normalize(propagation_array[0:4, 0:4]))

            propagated = propagation_array * masked
#            print(propagated[0:4, 0:4])
#            showimage(normalize(propagated[0:4, 0:4]))

#            showimage(equalize(propagated))
#            showimage(logscale(np.angle(propagated)))
#            showimage(equalize(propagated.real))

#            plt.cla()
#            plt.plot(np.sin((np.arange(50) / 3.)))
#            graph = fig2image(figure)
#            showimage(graph)

            reconstructed = get_ifft(propagated)
            intensity = get_intensity(reconstructed)
#            print(reconstructed[0:4, 0:4])
#            showimage(equalize(np.angle(reconstructed)))
#            showimage(equalize(intensity))
            showimage(equalize(get_shiftedfft(intensity)))
#            showimage(normalize(reconstructed.real))

#            showimage(normalize(masked[0:8, 0:8].real), 
#                normalize(propagation_array[0:8, 0:8].real),
#                normalize(reconstructed[0:8, 0:8].real),
#            )
        
#            showimage(normalize(masked[0:8, 0:8].imag), 
#                normalize(propagation_array[0:8, 0:8].imag),
#                normalize(reconstructed[0:8, 0:8].imag),
#            )
        
    return 0


if __name__ == "__main__":
    exit(main())
