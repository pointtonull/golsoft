#!/usr/bin/env python
#-*- coding: UTF-8 -*-

"""
Simple implementation of the Angular Spectrum Method to reconstruct lensless
holograms
"""

from automask import get_circles, get_holed_window, get_mask
from fmt import get_shiftedfft, get_ifft, get_shiftedifft
from image import equalize, get_intensity, get_centered, normalize
from numpy import exp, cos, sqrt
from scipy import optimize
from autopipe import showimage
import cache
import numpy as np


tau = 6.283185307179586
LAMBDA = 6.328e-07 # wave length
DX = 8.39e-6
DY = 8.46e-6
K = tau / LAMBDA # wave number
MASK_SOFTNESS = 2
MASK_R_SCALE = 3


@cache.hybrid
def apply_mask(array):
    """
    Try to filter out spurious data.
    """
    array = get_centered(array)
    shape = array.shape
    intensity = equalize(array)
    showimage(intensity)

    windowmaker = lambda x: np.kaiser(x, MASK_SOFTNESS)
    circles = sorted(get_circles(intensity, 3))
    print(circles)
    virtual_order, real_order, zero_order = circles

    centered = get_centered(array, real_order[1])

    window = get_holed_window(windowmaker, real_order[2] * MASK_R_SCALE, 10)
    mask = get_mask(shape, window)

    masked = get_centered(mask * centered)
    showimage(normalize(masked))
    return masked


#@cache.hybrid
def get_ref_beam(shape, alpha=90, beta=90):
    """
    Generate a reference beam array given the shape of the hologram and the
    directors angles
    """
    maxrow = shape[0] / 2
    maxcol = shape[1] / 2
    minrow, mincol = -maxrow, -maxcol
    row, col = np.ogrid[minrow:maxrow:1., mincol:maxcol:1.]
    cosa = cos(alpha)
    cosb = cos(beta)
    ref_beam = exp(1j * K * (cosa * col * DX + cosb * row * DY))
#    ref_beam = np.ones(shape)
    return ref_beam


#@cache.hybrid
def get_pea(hologram, distance, alpha=90, beta=90):
    """
    1. hologram x ref_beam
    2. shifted_fft(1)
    3. automask(2)
    4. center(3)
    5. propagation_factor_array(M) x 5
    6. shifted_ifft(5)
    """

    shape = hologram.shape
    ref_beam = get_ref_beam(shape, alpha, beta)
    rhologram = ref_beam * hologram

    frh = get_shiftedfft(rhologram)
    frh = get_centered(frh)
    masked = apply_mask(frh)

    maxrow = shape[0] / 2
    maxcol = shape[1] / 2
    minrow, mincol = -maxrow, -maxcol
    row, col = np.ogrid[minrow:maxrow:1., mincol:maxcol:1.]
    phase_correction_factor = sqrt(K * 1 - (LAMBDA * 232.7920143 * row) -
        (LAMBDA * 230.8658393 * col))
    propagation_array = exp(1j * phase_correction_factor * distance)
    propagation_array = get_centered(propagation_array)
    print("Propagation array")
    showimage(equalize(propagation_array.real))
    propagated = propagation_array * masked

    reconstructed = get_ifft(propagated)
    return reconstructed
    wrapped_phase = np.angle(reconstructed)
    return wrapped_phase


def get_strips_angle_radius(hologram):
    """
    Calculate the inclination angle and raiun for the given hologram strips
    """
    shape = hologram.shape
    diagonal = sum((dim ** 2 for dim in shape)) ** .5
    center = [dim / 2. for dim in shape]

    fft = get_intensity(get_shiftedfft(hologram))
    peak = get_circles(fft, 2)[1][1]
    peak = peak[0] - center[0], peak[1] - center[1]
    radius = (peak[0] ** 2 + peak[1] ** 2) ** .5 / diagonal * 2.
    angle = (1.570796325 - np.arctan2(*peak)) % 3.1415926536
    return angle, radius


#@cache.hybrid
def guess_angles(hologram):
    """
    Uses ML algoritms to guess the corrects directors angles for the given
    hologram
    """
    angle, radius = get_strips_angle_radius(hologram)
    ref_shape = [dim / 2 for dim in hologram.shape]

    def get_angles_fitness(args):
        """
        Learning fitness function
        """
        ref_beam = get_ref_beam(ref_shape, args[0], args[1])
        new_angle, new_radius = get_strips_angle_radius(ref_beam)
        angle_diff = (new_angle - angle) / 3.14159265
        radius_diff = new_radius - radius
        distance = (angle_diff ** 2 + radius_diff ** 2) ** .5
        return distance

    xinit = np.array([tau/4., tau/4.])
    optimizer = optimize.fmin # 66

    xend = optimizer(get_angles_fitness, xinit)
    return xend[0], xend[1]
