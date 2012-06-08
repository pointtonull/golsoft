#!/usr/bin/env python
#-*- coding: UTF-8 -*-

"""
Simple implementation of the Angular Spectrum Method to reconstruct lensless
holograms
"""

from autopipe import showimage
from automask import get_circles, get_holed_window, get_mask
from dft import get_shifted_dft, get_idft
from image import equalize, get_intensity, get_centered, normalize, logscale
from numpy import exp, cos, sqrt, sin
from scipy import optimize
import cache
import numpy as np


tau = 6.283185307179586
LAMBDA = 6.328e-07 # wave length
DX = 8.39e-6
DY = 8.46e-6
K = tau / LAMBDA # wave number
EPSILON = 1e-16


def get_auto_mask(spectrum, softness=0, radious_scale=1, zero_scale=1,
        cuttop=0):
    """
    Try to filter spurious data out.
    """
    shape = spectrum.shape
    intensity = get_intensity(spectrum)

    circles = sorted(get_circles(intensity, 3, 50))
    virtual_order, real_order, zero_order = circles
    peak_height, peak_center, peak_radious = real_order

    peak_radious = min([(abs(shape[0] / 3.5 - peak[2]), peak[2])
        for peak in circles])[1]

    windowmaker = lambda x: np.kaiser(x, softness)
    window = get_holed_window(windowmaker, peak_radious * radious_scale)
    mask = get_mask(shape, window, peak_center)

    zerowindow = get_holed_window(windowmaker, peak_radious * zero_scale)
    zeromask = 1 - get_mask(shape, zerowindow, zero_order[1])
    mask *= zeromask

    masked_intensity = mask * intensity

    cutoff = masked_intensity > (masked_intensity.max()
        - masked_intensity.ptp() * cuttop)
    mask[cutoff] = 0
    masked = mask * spectrum

    centered = get_centered(intensity, peak_center)
    masked = get_centered(masked, peak_center)
    mask = get_centered(mask, peak_center)

    return mask, masked, centered


def get_ref_beam(shape, cos_alpha=EPSILON, cos_beta=EPSILON):
    """
    Generate a reference beam array given the shape of the hologram and the
    directors angles
    """
    maxrow = shape[0] / 2
    maxcol = shape[1] / 2
    minrow, mincol = -maxrow, -maxcol
    row, col = np.ogrid[minrow:maxrow:1., mincol:maxcol:1.]

    ref_beam = exp(1j * K * (cos_alpha * col * DX + cos_beta * row * DY))
    
    return ref_beam


def get_pea(hologram, distance, cos_alpha=EPSILON, cos_beta=EPSILON,
        radious_scale=1, softness=1):
    """
    1. hologram x ref_beam
    2. shifted_fft(1)
    3. automask(2)
    4. center(3)
    5. propagation_factor_array(M) x 5
    6. shifted_ifft(5)
    """

    shape = hologram.shape
    ref_beam = get_ref_beam(shape, cos_alpha, cos_beta)
    rhologram = ref_beam * hologram

    frh = get_shifted_dft(rhologram)
    mask, masked, centered = get_auto_mask(get_intensity(frh), softness, 
        radious_scale)
    masked = frh * mask

    maxrow = shape[0] / 2
    maxcol = shape[1] / 2
    minrow, mincol = -maxrow, -maxcol
    row, col = np.ogrid[minrow:maxrow:1., mincol:maxcol:1.]
    propagation_array = get_propagation_array(shape, distance)
    propagated = propagation_array * masked

    reconstructed = get_idft(propagated)
    return reconstructed
 

@cache.hybrid(reset=0)
def get_propagation_array(shape, distance):
    maxrow = shape[0] / 2.
    maxcol = shape[1] / 2.
    minrow, mincol = -maxrow, -maxcol
    row, col = np.ogrid[minrow:maxrow:1, mincol:maxcol:1]
    phase_correction_factor = K * sqrt(1 - (LAMBDA * 232.7920143 * row) ** 2
        - (LAMBDA * 230.8658393 * col) ** 2)
    propagation_array = exp(1j * phase_correction_factor * distance)
    return propagation_array


def get_distance(point1, point2):
    distance = ((point1[0] - point2[0]) ** 2
        + (point1[1] - point2[1]) ** 2) ** .5
    return distance


def get_peak_coords(hologram):
    """
    """
    shape = hologram.shape
    center = [dim / 2. for dim in shape]

    dft = get_intensity(get_shifted_dft(hologram.real))
    circles = [(-get_distance(center, circle[1]), circle[0], circle[1])
        for circle in get_circles(dft, 2, 20)]
    circles.sort()
    peak = circles[0][2]
    peaks_row = (peak[0] - center[0]) / float(center[0])
    peaks_col = (peak[1] - center[1]) / float(center[1])
    return peaks_row, peaks_col


@cache.hybrid(reset=0)
def get_refbeam_peak_coords(alpha, beta):
    ref_beam = get_ref_beam((256, 256), alpha, beta)
    row, col = get_peak_coords(ref_beam)
    return row, col


@cache.hybrid(reset=0)
def guess_director_cosines(hologram):
    """
    guess the optimums directors angles for the given hologram
    """
    ref_peak = get_peak_coords(hologram)
    print("ref_peak: %4.3f, %4.3f" % ref_peak)

    def get_fitness(args):
        """
        Learning fitness function
        """
        new_peak = get_refbeam_peak_coords(*args)
        distance = get_distance(ref_peak, new_peak)
        return distance

    xinit = np.array([tau / 4, tau /4])
    xend = generic_minimizer(get_fitness, xinit)
    return xend[0], xend[1]


def generic_minimizer(fitness_func, initial_guess, optimizers=None):
    """
    A common interface to various minimization algorithms
    """

    if optimizers == None:
        optimizers = [
            optimize.fmin, # 66
            optimize.fmin_powell,
        ]

    best_result = None
    for optimizer in optimizers:
        xend = optimizer(fitness_func, initial_guess, disp=False)
        last_result = fitness_func(xend)
        if best_result is None or last_result < best_result:
            best_guess = xend
            best_result = last_result

    return best_guess


@cache.hybrid(reset=0)
def calculate_director_cosines(hologram):
    """
    Calculate the director cosines using the spectral proyection formula
    """
    peak = get_peak_coords(hologram)
    freq_rows, freq_cols = peak
    freq_rows /= 2 * DY
    freq_cols /= 2 * DX
    cos_alpha = freq_cols * LAMBDA
    cos_beta = freq_rows * LAMBDA

    return cos_alpha, cos_beta
