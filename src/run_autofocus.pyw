<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> experimental
#!/usr/bin/env python
#-*- coding: UTF-8 -*-

import sys

from scipy import misc, optimize, stats
import matplotlib.pyplot as plt
import numpy as np

from automask import get_mask
from autopipe import showimage
<<<<<<< HEAD
from dft import get_shifted_dft, get_idft
from image import imread, normalize, get_intensity, equalize
from pea import get_auto_mask, generic_minimizer
from pea import get_propagation_array, get_distance
from ranges import frange
import cache

=======
from dft import get_shifted_dft, get_shifted_idft
from image import imread, normalize, get_intensity, equalize
from pea import get_auto_mask, generic_minimizer, angle2
from pea import get_propagation_array, get_distance
from ranges import frange
import cache
#from fresnel import get_chirp, get_wrapped_formula_array
>>>>>>> experimental


class Methods(list):
    def __call__(self, func):
        self.append(func)
        return func
methods = Methods()


def get_highpass_mask(shape, radius=0.2, softness=0):
    radius = round(min(shape) * (1 - radius))
    window = np.kaiser(radius, softness)
    mask = 1 - get_mask(shape, window)
    return mask


def get_lowpass_mask(shape, radius=0.2, softness=0):
    radius = round(min(shape) * radius)
    window = np.kaiser(radius, softness)
    mask = get_mask(shape, window)
    return mask

<<<<<<< HEAD

=======
>>>>>>> experimental
@methods
@cache.hybrid
def get_var(masked_spectrum, distance):
    propagation_array = get_propagation_array(masked_spectrum.shape, distance)
    propagated = propagation_array * masked_spectrum
<<<<<<< HEAD
    reconstructed = get_idft(propagated)
#    intensity = get_intensity(reconstructed)
    intensity = np.abs(reconstructed)
    fitness = intensity.var()
=======
    reconstructed = get_shifted_idft(propagated)
    module = np.abs(reconstructed)
    fitness = module.var()
>>>>>>> experimental
    return fitness

@methods
def get_diff_var(masked_spectrum, distance):
    fitness = get_var(masked_spectrum, distance)
    fitness -= get_var(masked_spectrum, -distance)
    return fitness


<<<<<<< HEAD
=======
@methods
@cache.hybrid
def get_int_var(masked_spectrum, distance):
    propagation_array = get_propagation_array(masked_spectrum.shape, distance)
    propagated = propagation_array * masked_spectrum
    reconstructed = get_shifted_idft(propagated)
    intensity = get_intensity(reconstructed)
    fitness = intensity.var()
    return fitness

@methods
def get_diff_int_var(masked_spectrum, distance):
    fitness = get_int_var(masked_spectrum, distance)
    fitness -= get_int_var(masked_spectrum, -distance)
    return fitness


>>>>>>> experimental
#@methods
@cache.hybrid
def get_lowpass_var(masked_spectrum, distance):
    propagation_array = get_propagation_array(masked_spectrum.shape, distance)
    propagated = propagation_array * masked_spectrum
    lowpass_mask = get_lowpass_mask(propagated.shape, .4)
    propagated = lowpass_mask * propagated
<<<<<<< HEAD
    reconstructed = get_idft(propagated)
=======
    reconstructed = get_shifted_idft(propagated)
>>>>>>> experimental
    intensity = get_intensity(reconstructed)
    fitness = intensity.var()
    return fitness

#@methods
def get_diff_lowpass_var(masked_spectrum, distance):
    fitness = get_lowpass_var(masked_spectrum, distance)
    fitness -= get_lowpass_var(masked_spectrum, -distance)
    return fitness



<<<<<<< HEAD
<<<<<<< HEAD
#@methods
def phase_detection(masked_spectrum, distance):
    shape = masked_spectrum.shape
    shape_center = [dim / 2 for dim in shape]

    filtered_hologram = get_shifted_idft(masked_spectrum)
    focus_mask = get_mask(shape, np.ones(20), shape_center)
    focus_feature = filtered_hologram * focus_mask
    feature_spectrum = get_shifted_dft(focus_feature)

    propagation_array = get_propagation_array(shape, distance)
    propagated_spectrum = propagation_array * feature_spectrum

    propagated_hologram = get_shifted_idft(propagated_spectrum)

    radious = 50
    separation = 400

    window = np.ones(radious)

    spectrum2sensor = np.angle
    spectrum2sensor = np.abs
    spectrum2sensor = get_intensity
    spectrum2sensor = normalize

    left_center = shape_center[0], shape_center[1] - separation / 2
    left_mask = get_mask(shape, window, left_center)
#    left_bundle = left_mask * propagated_spectrum
    left_bundle = left_mask * propagated_hologram
    left_sensor = spectrum2sensor(get_shifted_dft(left_bundle))
    left_peak = get_circles(left_sensor, 1, 50)[0][1]

    right_center = shape_center[0], shape_center[1] + separation / 2
    right_mask = get_mask(shape, window, right_center)
    right_bundle = right_mask * propagated_hologram
    right_sensor = spectrum2sensor(get_shifted_dft(right_bundle))
    right_peak = get_circles(right_sensor, 1, 50)[0][1]

#    showimage(normalize(np.maximum(left_bundle, right_bundle)))
    showimage(normalize(np.maximum(left_sensor, right_sensor)))

    distance = get_distance(left_peak, right_peak)
    fitness = distance
    return fitness


=======
>>>>>>> experimental
=======
>>>>>>> experimental
#@methods
@cache.hybrid(reset=0)
def get_highpass_var(masked_spectrum, distance):
    propagation_array = get_propagation_array(masked_spectrum.shape, distance)
    propagated = propagation_array * masked_spectrum
    highpass_mask = get_highpass_mask(propagated.shape, .4)
    propagated = highpass_mask * propagated
<<<<<<< HEAD
    reconstructed = get_idft(propagated)
=======
    reconstructed = get_shifted_idft(propagated)
>>>>>>> experimental
    intensity = get_intensity(reconstructed)
    fitness = intensity.var()
    return fitness

#@methods
def get_diff_highpass_var(masked_spectrum, distance):
    fitness = get_highpass_var(masked_spectrum, distance)
    fitness -= get_highpass_var(masked_spectrum, -distance)
    return fitness



#@methods
@cache.hybrid
def get_var_over_hpass_var(*args):
    fitness = get_var(*args) / get_highpass_var(*args)
    return fitness

#@methods
def get_diff_var_over_hpass_var(masked_spectrum, distance):
    fitness = get_var_over_hpass_var(masked_spectrum, distance)
    fitness -= get_var_over_hpass_var(masked_spectrum, -distance)
    return fitness




#@methods
@cache.hybrid
def get_lpass_var_over_hpass_var(*args):
    fitness = get_lowpass_var(*args) / get_highpass_var(*args)
    return fitness


#@methods
def get_diff_lpass_var_over_hpass_var(masked_spectrum, distance):
    fitness = get_lpass_var_over_hpass_var(masked_spectrum, distance)
    fitness -= get_lpass_var_over_hpass_var(masked_spectrum, -distance)
    return fitness


<<<<<<< HEAD
def get_best_contrast_zone(hologram, shape=(400, 400)):
=======
def get_best_contrast_zone(hologram, shape=(256, 256)):
>>>>>>> experimental
    assert shape[0] <= hologram.shape[0]
    assert shape[1] <= hologram.shape[1]
    rows = hologram.shape[0] - shape[0] + 1
    cols = hologram.shape[1] - shape[1] + 1
    rowsvar = hologram.var(0)
    colsvar = hologram.var(1)
    sumsrowsranges = np.array([rowsvar[top:top + shape[0]].sum()
        for top in xrange(rows)])
    sumscolsranges = np.array([colsvar[left:left + shape[1]].sum()
        for left in xrange(cols)])
    toprow = sumsrowsranges.argmax()
    leftcol = sumscolsranges.argmax()
    print(rows, cols)
    return hologram[toprow:toprow + shape[0], leftcol:leftcol + shape[1]]


def guess_focus_distance(masked_spectrum, extractor):

    def fitness(args):
        return extractor(masked_spectrum, args)

    results = []
    for distance in frange(0, .15, 3):
        xend = generic_minimizer(fitness, distance, [optimize.fmin])
        results.append(xend)
    return results



def main():
    images = [(filename, imread(filename, True))
        for filename in sys.argv[1:]]
    if len(images) < 2:
        if not images:
            lena = misc.imresize(misc.lena(), .5)
            images = [("lena", lena)]

    figure = plt.figure()

<<<<<<< HEAD
    graph_distances = [distance for distance in frange(0.0, 2**-2, 80)]
=======
    graph_distances = [distance for distance in frange(0.0, 2**-2, 160)]
>>>>>>> experimental

    for filename, hologram in images:
        print("\nOriginal image: %s" % filename)
        shape = hologram.shape
        showimage(hologram)

        best_zone = get_best_contrast_zone(hologram)
        showimage(best_zone)
        print(best_zone.shape)

        spectrum = get_shifted_dft(hologram)
        mask, masked_spectrum, centered = get_auto_mask(spectrum,
            softness=0, radious_scale=1.5)
        showimage(equalize(centered), equalize(masked_spectrum))

        zone_spectrum = get_shifted_dft(best_zone)
<<<<<<< HEAD
        zone_spectrum = spectrum
=======
>>>>>>> experimental
        mask, zone_masked_spectrum, centered = get_auto_mask(zone_spectrum,
            softness=0, radious_scale=1.0)
        showimage(equalize(centered), equalize(zone_masked_spectrum))

        for method in methods:
            print("\nMethod: %s\n" % method.func_name)

            fitness_values = [method(zone_masked_spectrum, distance)
                for distance in graph_distances]
            plt.cla()
            plt.plot(graph_distances, fitness_values, c="blue")

            localmins = guess_focus_distance(zone_masked_spectrum, method)
            fitness = [method(zone_masked_spectrum, dst)
                for dst in localmins]

            plt.scatter(localmins, fitness, c="green")

            bestfitness, globalmin = min(zip(fitness, localmins))
            plt.scatter(globalmin, bestfitness, c="red")

            showimage(figure)
<<<<<<< HEAD

            propagation_array = get_propagation_array(shape, globalmin)
            propagated = masked_spectrum * propagation_array
            reconstructed = get_idft(propagated)
            showimage(normalize(np.abs(reconstructed)), 
                normalize(np.arctan2(reconstructed.real,
                reconstructed.imag)))
=======
            
            propagation_array = get_propagation_array(zone_spectrum.shape, globalmin)
            propagated = zone_masked_spectrum * propagation_array
            reconstructed = get_shifted_idft(propagated)
            showimage(normalize(np.abs(reconstructed)), 
                normalize(angle2(reconstructed)))
            print(propagated.shape, reconstructed.shape)

            propagation_array = get_propagation_array(shape, globalmin)
            propagated = masked_spectrum * propagation_array
            reconstructed = get_shifted_idft(propagated)
            showimage(normalize(np.abs(reconstructed)), 
                normalize(angle2(reconstructed)))
>>>>>>> experimental
            print globalmin, "\n"

    return 0


if __name__ == "__main__":
    exit(main())
<<<<<<< HEAD
=======
#!/usr/bin/env python
#-*- coding: UTF-8 -*-

import sys

from scipy import misc, optimize, stats
import matplotlib.pyplot as plt
import numpy as np

from automask import get_mask
from autopipe import showimage
from dft import get_shifted_dft, get_shifted_idft
from image import imread, normalize, get_intensity, equalize
from pea import get_auto_mask, generic_minimizer, angle2
from pea import get_propagation_array, get_distance
from ranges import frange
import cache
#from fresnel import get_chirp, get_wrapped_formula_array


class Methods(list):
    def __call__(self, func):
        self.append(func)
        return func
methods = Methods()


def get_highpass_mask(shape, radius=0.2, softness=0):
    radius = round(min(shape) * (1 - radius))
    window = np.kaiser(radius, softness)
    mask = 1 - get_mask(shape, window)
    return mask


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
    reconstructed = get_shifted_idft(propagated)
    module = np.abs(reconstructed)
    fitness = module.var()
    return fitness

@methods
def get_diff_var(masked_spectrum, distance):
    fitness = get_var(masked_spectrum, distance)
    fitness -= get_var(masked_spectrum, -distance)
    return fitness


@methods
@cache.hybrid
def get_int_var(masked_spectrum, distance):
    propagation_array = get_propagation_array(masked_spectrum.shape, distance)
    propagated = propagation_array * masked_spectrum
    reconstructed = get_shifted_idft(propagated)
    intensity = get_intensity(reconstructed)
    fitness = intensity.var()
    return fitness

@methods
def get_diff_int_var(masked_spectrum, distance):
    fitness = get_int_var(masked_spectrum, distance)
    fitness -= get_int_var(masked_spectrum, -distance)
    return fitness


#@methods
@cache.hybrid
def get_lowpass_var(masked_spectrum, distance):
    propagation_array = get_propagation_array(masked_spectrum.shape, distance)
    propagated = propagation_array * masked_spectrum
    lowpass_mask = get_lowpass_mask(propagated.shape, .4)
    propagated = lowpass_mask * propagated
    reconstructed = get_shifted_idft(propagated)
    intensity = get_intensity(reconstructed)
    fitness = intensity.var()
    return fitness

#@methods
def get_diff_lowpass_var(masked_spectrum, distance):
    fitness = get_lowpass_var(masked_spectrum, distance)
    fitness -= get_lowpass_var(masked_spectrum, -distance)
    return fitness



#@methods
@cache.hybrid(reset=0)
def get_highpass_var(masked_spectrum, distance):
    propagation_array = get_propagation_array(masked_spectrum.shape, distance)
    propagated = propagation_array * masked_spectrum
    highpass_mask = get_highpass_mask(propagated.shape, .4)
    propagated = highpass_mask * propagated
    reconstructed = get_shifted_idft(propagated)
    intensity = get_intensity(reconstructed)
    fitness = intensity.var()
    return fitness

#@methods
def get_diff_highpass_var(masked_spectrum, distance):
    fitness = get_highpass_var(masked_spectrum, distance)
    fitness -= get_highpass_var(masked_spectrum, -distance)
    return fitness



#@methods
@cache.hybrid
def get_var_over_hpass_var(*args):
    fitness = get_var(*args) / get_highpass_var(*args)
    return fitness

#@methods
def get_diff_var_over_hpass_var(masked_spectrum, distance):
    fitness = get_var_over_hpass_var(masked_spectrum, distance)
    fitness -= get_var_over_hpass_var(masked_spectrum, -distance)
    return fitness




#@methods
@cache.hybrid
def get_lpass_var_over_hpass_var(*args):
    fitness = get_lowpass_var(*args) / get_highpass_var(*args)
    return fitness


#@methods
def get_diff_lpass_var_over_hpass_var(masked_spectrum, distance):
    fitness = get_lpass_var_over_hpass_var(masked_spectrum, distance)
    fitness -= get_lpass_var_over_hpass_var(masked_spectrum, -distance)
    return fitness


def get_best_contrast_zone(hologram, shape=(256, 256)):
    assert shape[0] <= hologram.shape[0]
    assert shape[1] <= hologram.shape[1]
    rows = hologram.shape[0] - shape[0] + 1
    cols = hologram.shape[1] - shape[1] + 1
    rowsvar = hologram.var(0)
    colsvar = hologram.var(1)
    sumsrowsranges = np.array([rowsvar[top:top + shape[0]].sum()
        for top in xrange(rows)])
    sumscolsranges = np.array([colsvar[left:left + shape[1]].sum()
        for left in xrange(cols)])
    toprow = sumsrowsranges.argmax()
    leftcol = sumscolsranges.argmax()
    print(rows, cols)
    return hologram[toprow:toprow + shape[0], leftcol:leftcol + shape[1]]


def guess_focus_distance(masked_spectrum, extractor):

    def fitness(args):
        return extractor(masked_spectrum, args)

    results = []
    for distance in frange(0, .15, 3):
        xend = generic_minimizer(fitness, distance, [optimize.fmin])
        results.append(xend)
    return results



def main():
    images = [(filename, imread(filename, True))
        for filename in sys.argv[1:]]
    if len(images) < 2:
        if not images:
            lena = misc.imresize(misc.lena(), .5)
            images = [("lena", lena)]

    figure = plt.figure()

    graph_distances = [distance for distance in frange(0.0, 2**-2, 160)]

    for filename, hologram in images:
        print("\nOriginal image: %s" % filename)
        shape = hologram.shape
        showimage(hologram)

        best_zone = get_best_contrast_zone(hologram)
        showimage(best_zone)
        print(best_zone.shape)

        spectrum = get_shifted_dft(hologram)
        mask, masked_spectrum, centered = get_auto_mask(spectrum,
            softness=0, radious_scale=1.5)
        showimage(equalize(centered), equalize(masked_spectrum))

        zone_spectrum = get_shifted_dft(best_zone)
        mask, zone_masked_spectrum, centered = get_auto_mask(zone_spectrum,
            softness=0, radious_scale=1.0)
        showimage(equalize(centered), equalize(zone_masked_spectrum))

        for method in methods:
            print("\nMethod: %s\n" % method.func_name)

            fitness_values = [method(zone_masked_spectrum, distance)
                for distance in graph_distances]
            plt.cla()
            plt.plot(graph_distances, fitness_values, c="blue")

            localmins = guess_focus_distance(zone_masked_spectrum, method)
            fitness = [method(zone_masked_spectrum, dst)
                for dst in localmins]

            plt.scatter(localmins, fitness, c="green")

            bestfitness, globalmin = min(zip(fitness, localmins))
            plt.scatter(globalmin, bestfitness, c="red")

            showimage(figure)
            
            propagation_array = get_propagation_array(zone_spectrum.shape, globalmin)
            propagated = zone_masked_spectrum * propagation_array
            reconstructed = get_shifted_idft(propagated)
            showimage(normalize(np.abs(reconstructed)), 
                normalize(angle2(reconstructed)))
            print(propagated.shape, reconstructed.shape)

            propagation_array = get_propagation_array(shape, globalmin)
            propagated = masked_spectrum * propagation_array
            reconstructed = get_shifted_idft(propagated)
            showimage(normalize(np.abs(reconstructed)), 
                normalize(angle2(reconstructed)))
            print globalmin, "\n"

    return 0


if __name__ == "__main__":
    exit(main())
>>>>>>> experimental
=======
>>>>>>> experimental
