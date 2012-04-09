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
from pea import apply_mask, generic_minimizer
from pea import guess_director_angles, get_ref_beam, get_propagation_array
from scipy import misc, ndimage
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


def get_fitness(masked_spectrum, distance):
    """
    Learning fitness function
    """
    propagation_array = get_propagation_array(masked_spectrum.shape, distance)
    propagated = propagation_array * masked_spectrum
    reconstructed = get_ifft(propagated)
    intensity = get_intensity(reconstructed)
#    showimage(equalize(intensity))
    diff = np.diff(intensity)
    fitness = ndimage.variance(intensity)
    print(distance, fitness)
    return fitness


def guess_focus_distance(masked_spectrum):
    shape = masked_spectrum.shape

    get_fitness = lambda args: get_ptp(masked_spectrum, args[0])

    xinit = np.array([0])
    xend = generic_minimizer(get_fitness, xinit)
    return xend


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
        masked_spectrum = apply_mask(spectrum, softness=0, radius_scale=3)

#        for distance in frange(-0.0359375 , 2**-4, 5):

#        if True:
#            distance = guess_focus_distance(masked_spectrum)

#            print(distance)
#            print(get_ptp(masked_spectrum, distance))

#            propagation_array = get_propagation_array(shape, distance)
#            propagated = propagation_array * masked_spectrum
#            showimage(equalize(propagated.imag))

#            showimage(equalize(propagated))
#            showimage(logscale(np.angle(propagated)))
#            showimage(equalize(propagated.real))

        distances = [distance for distance in frange(-0.05 , 2**-2, 100)]
        fitness_values = [get_fitness(masked_spectrum, distance)
            for distance in distances]
        plt.cla()
#        plt.plot(distances, ptp_values)
        plt.scatter(distances, fitness_values)
        graph = fig2image(figure)
        showimage(graph)

#            reconstructed = get_ifft(propagated)

#            showimage(equalize(np.angle(reconstructed))) # phase

#            intensity = get_intensity(reconstructed)
#            showimage(equalize(intensity)) # module

#            showimage(equalize(get_shiftedfft(intensity)))

    return 0


if __name__ == "__main__":
    exit(main())
