#!/usr/bin/env python
#-*- coding: UTF-8 -*-

import sys

from scipy import misc, optimize
import matplotlib.pyplot as plt
import numpy as np

from automask import get_mask
from autopipe import showimage
from dft import get_shifted_dft, get_shifted_idft
from image import imread, normalize, get_intensity, equalize
from pea import get_auto_mask, generic_minimizer, angle2
from pea import get_propagation_array
from ranges import frange
import cache


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


@cache.hybrid
def guess_focus_distance(masked_spectrum, extractor=get_diff_var):

    def fitness(args):
        return extractor(masked_spectrum, args)

    results = []
    for distance in frange(0, .15, 3):
        xend = generic_minimizer(fitness, distance, [optimize.fmin])
        results.append((fitness(xend), xend[0]))
    return sorted(results)[0][1]
