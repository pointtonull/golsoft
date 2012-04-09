#!/usr/bin/env python
#-*- coding: UTF-8 -*-

"""
Simple implementation of the Angular Spectrum Method to reconstruct lensless
holograms
"""

#from autopipe import showimage
from automask import get_circles, get_holed_window, get_mask
from fmt import get_shiftedfft, get_ifft, get_shiftedifft
from image import equalize, get_intensity, get_centered, normalize, logscale
from numpy import exp, cos, sqrt
from scipy import optimize
import cache
import numpy as np


tau = 6.283185307179586
LAMBDA = 6.328e-07 # wave length
DX = 8.39e-6
DY = 8.46e-6
K = tau / LAMBDA # wave number
MASK_SOFTNESS = 2
MASK_R_SCALE = 2 # FIXME


#@cache.hybrid
def apply_mask(array):
    """
    Try to filter out spurious data.
    """
    array = get_centered(array)
    shape = array.shape
    intensity = equalize(array)
#    showimage(intensity)

    windowmaker = lambda x: np.kaiser(x, 0)#FIXME
    circles = sorted(get_circles(intensity, 3, 50))
    print(circles)
    virtual_order, real_order, zero_order = circles

    centered = get_centered(array, real_order[1])
#    showimage(equalize(centered))

    window = get_holed_window(windowmaker, real_order[2] * MASK_R_SCALE,
        0)#FIXME
    mask = get_mask(shape, window)
#    showimage(normalize(mask))

    masked = get_centered(mask * centered)
#    showimage(logscale(masked))
    return masked


def get_ref_beam(shape, alpha=tau/4, beta=tau/4):
    """
    Generate a reference beam array given the shape of the hologram and the
    directors angles
    """
    cos_alpha = cos(alpha)
    cos_beta = cos(beta)
    maxrow = shape[0] / 2
    maxcol = shape[1] / 2
    minrow, mincol = -maxrow, -maxcol
    row, col = np.ogrid[minrow:maxrow:1., mincol:maxcol:1.]
    ref_beam = exp(1j * K * (cos_alpha * col * DX + cos_beta * row * DY))
    return ref_beam


#@cache.hybrid
def get_pea(hologram, distance, alpha=tau/4, cos_beta=tau/4):
    """
    1. hologram x ref_beam
    2. shifted_fft(1)
    3. automask(2)
    4. center(3)
    5. propagation_factor_array(M) x 5
    6. shifted_ifft(5)
    """

    cos_alpha = cos(alpha)
    cos_beta = cos(beta)
    shape = hologram.shape
    ref_beam = get_ref_beam(shape, alpha, beta)
    rhologram = ref_beam * hologram

    frh = get_shiftedfft(rhologram)
    masked = apply_mask(frh)

    maxrow = shape[0] / 2
    maxcol = shape[1] / 2
    minrow, mincol = -maxrow, -maxcol
    row, col = np.ogrid[minrow:maxrow:1., mincol:maxcol:1.]
    phase_correction_factor = K * sqrt(1 - (LAMBDA * 232.7920143 * row)**2 -
        (LAMBDA * 230.8658393 * col)**2)
    propagation_array = exp(1j * phase_correction_factor * distance)
    print("Propagation array")
#    showimage(equalize(propagation_array.real))
    propagated = propagation_array * masked

    reconstructed = get_ifft(propagated)
    return reconstructed
    wrapped_phase = np.angle(reconstructed)
    return wrapped_phase


def get_distance(point1, point2):
    distance = ((point1[0] - point2[0]) ** 2
        + (point1[1] - point2[1]) ** 2) ** .5
    return distance


def get_peak_coords(hologram):
    """
    """
    shape = hologram.shape
    center = [dim / 2. for dim in shape]

    fft = get_intensity(get_shiftedfft(hologram.real))
    circles = [(-get_distance(center, circle[1]), circle[0], circle[1])
        for circle in get_circles(fft, 2, 20)]
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
def guess_director_angles(hologram):
    """
    guess the optimums directors angles for the given hologram
    """
    ref_peak = get_peak_coords(hologram)

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


def generic_minimizer(fitness_func, initial_guess, epsilon=5e-3):
    optimizers = [
        optimize.fmin, # 66
        optimize.fmin_powell,
#        optimize.fmin_cg,
#        optimize.fmin_slsqp,
#        optimize.fmin_bfgs,
    ]

    best_result = None
    for optimizer in optimizers:
        xend = optimizer(fitness_func, initial_guess, disp=False)
        last_result = fitness_func(xend)
        if best_result is None or last_result < best_result:
            best_guess = xend
            best_result = last_result
            print(optimizer.func_name, best_guess, last_result)
            if last_result < epsilon:
                break

    return best_guess
