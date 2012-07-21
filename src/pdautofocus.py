#!/usr/bin/env python
#-*- coding: UTF-8 -*-

import numpy as np

from automask import get_mask, get_circles
from image import normalize, get_intensity
from dft import get_shifted_dft, get_shifted_idft
from pea import get_propagation_array, get_distance

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
#    showimage(normalize(np.maximum(left_sensor, right_sensor)))

    distance = get_distance(left_peak, right_peak)
    fitness = distance
    return fitness

