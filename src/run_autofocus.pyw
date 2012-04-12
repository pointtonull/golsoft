#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from autopipe import showimage
from dft import get_shifted_dft, get_idft
from image import imread, normalize, logscale, get_intensity
from pea import apply_mask, generic_minimizer
from pea import guess_director_angles, get_ref_beam, get_propagation_array
from ranges import frange
from scipy import misc, ndimage
#import Image
import matplotlib.pyplot as plt
import numpy as np
import sys



def get_fitness(masked_spectrum, distance):
    """
    Learning fitness function
    """
    propagation_array = get_propagation_array(masked_spectrum.shape, distance)
    propagated = propagation_array * masked_spectrum
    reconstructed = get_idft(propagated)
    intensity = get_intensity(reconstructed)
    fitness = ndimage.variance(intensity)
    return fitness


def guess_focus_distance(masked_spectrum):

    def fitness(args):
        return get_fitness(masked_spectrum, args)

    xinit = np.array([0])
    xend = generic_minimizer(fitness, xinit)
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
        image = normalize(image)
        showimage(image)

        alpha, beta = guess_director_angles(image)
        ref_beam = get_ref_beam(shape, alpha, beta)
        rhologram = ref_beam * image

        spectrum = get_shifted_dft(rhologram)
        masked_spectrum = apply_mask(spectrum, softness=0, radius_scale=3)

        distances = [distance for distance in frange(-0.05 , 2**-2, 100)]
        fitness_values = [get_fitness(masked_spectrum, distance)
            for distance in distances]
        plt.cla()
        plt.scatter(distances, fitness_values)

        distance = guess_focus_distance(masked_spectrum)
        if abs(distance) > 2:
            distance = 0

        fitness = get_fitness(masked_spectrum, distance)
        plt.scatter(distance, fitness, c="red")
        showimage(figure)
        print("Calculated distance: %1.5fm" % distance)

#        reconstructed = get_idft(propagated)
#        showimage(equalize(np.angle(reconstructed))) # phase
#        intensity = get_intensity(reconstructed)
#        showimage(equalize(intensity)) # module
#        showimage(equalize(get_shifted_dft(intensity)))

    return 0


if __name__ == "__main__":
    exit(main())
