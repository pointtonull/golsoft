#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from autopipe import showimage
from dft import get_shifted_dft, get_idft
from image import imread, normalize, get_intensity, equalize, get_centered
from pea import apply_mask, generic_minimizer
from pea import calculate_director_cosines, get_ref_beam, get_propagation_array
from fmt import get_mask
from ranges import frange
import cache
from scipy import misc, ndimage, optimize, stats
#import Image
import matplotlib.pyplot as plt
import numpy as np
import sys


class Methods(list):
    def __call__(self, func):
        self.append(func)
        return func

methods = Methods()


@cache.hybrid
def get_highpass_mask(shape, radius=0.2, softness=0):
    radius = round(min(shape) * (1 - radius))
    window = np.kaiser(radius, softness)
    mask = 1 - get_mask(shape, window)
    return mask


@cache.hybrid
def get_lowpass_mask(shape, radius=0.2, softness=0):
    radius = round(min(shape) * radius)
    window = np.kaiser(radius, softness)
    mask = get_mask(shape, window)
    return mask


@methods
@cache.hybrid
def get_var(masked_spectrum, distance):
    propagation_array = get_propagation_array(masked_spectrum.shape, distance)
    propagated = propagation_array * masked_spectrum
    reconstructed = get_idft(propagated)
    intensity = get_intensity(reconstructed)
    fitness = intensity.var()
    return fitness


#@methods
@cache.hybrid
def get_lowpass_var(masked_spectrum, distance):
    propagation_array = get_propagation_array(masked_spectrum.shape, distance)
    propagated = propagation_array * masked_spectrum
    lowpass_mask = get_lowpass_mask(propagated.shape, .4)
    propagated = lowpass_mask * propagated
    reconstructed = get_idft(propagated)
    intensity = get_intensity(reconstructed)
    fitness = intensity.var()
    return fitness



#@methods
@cache.hybrid(reset=0)
def get_highpass_var(masked_spectrum, distance):
    propagation_array = get_propagation_array(masked_spectrum.shape, distance)
    propagated = propagation_array * masked_spectrum
    highpass_mask = get_highpass_mask(propagated.shape, .4)
    propagated = highpass_mask * propagated
    reconstructed = get_idft(propagated)
    intensity = get_intensity(reconstructed)
    fitness = intensity.var()
    return fitness



@methods
@cache.hybrid
def get_var_over_hpass_var(*args):
    fitness = get_var(*args) / get_highpass_var(*args)
    return fitness



@methods
@cache.hybrid
def get_lpass_var_over_hpass_var(*args):
    fitness = get_lowpass_var(*args) / get_highpass_var(*args)
    return fitness


def get_groups(values, number):
    if number == 1 or len(values) == 0:
        yield [values]
    else:
        for cant in range(1, len(values) - number + 2):
            left = [values[:cant]]
            remainder = values[cant:]
            for right in get_groups(remainder, number - 1):
                yield left + right 


def get_metavar(groups):
    metavar = sum((np.var(group) for group in groups))
    return metavar


def autogroup(values):
    if len(values) > 2:
        values = sorted(values)
        groupings = [groups for cant in range(1, len(values) + 1)
            for groups in get_groups(values, cant)]
        lengths = [len(groups) for groups in groupings]
        variances = [get_metavar(groups) for groups in groupings]
        slope, intercept, r_value = stats.linregress(lengths, variances)[:3]
        regresion = [length * slope + intercept for length in lengths]
        diffs = np.array(variances) - np.array(regresion)
        pos = diffs.argmin()
        return groupings[pos]
    else:
        return values


def guess_focus_distance(masked_spectrum, extractor):

    def fitness(args):
        return extractor(masked_spectrum, args)

    results = []
    for distance in frange(0, .15, 6):
        xend = generic_minimizer(fitness, distance, [optimize.fmin])
        results.append(xend)
    return results


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
