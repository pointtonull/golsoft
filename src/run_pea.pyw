# !/usr/bin/env python
# -*- coding: UTF-8 -*-


from automask import get_circles, get_holed_window, get_mask
from autopipe import showimage, blue
import cache
from image import equalize, logscale, get_intensity, get_centered, imread
from fmt import get_shiftedfft, get_shiftedifft, get_fft, get_ifft
from numpy import exp, cos, sin, sqrt, angle
from random import sample
from scipy import misc, ndimage
import numpy as np
import sys


tau = 6.283185307179586
LAMBDA = 6.328e-07 # wave length
DX = 8.39e-6
DY = 8.46e-6
K = tau / LAMBDA # wave number
MASK_SOFTNESS = 2
MASK_R_SCALE = 3


def apply_mask(array):
    oarray = array
    array = get_centered(array) # TODO: validate (ACM)
    shape = array.shape
    intensity = equalize(array)

    print("Intensity map on apply_mask:")
    showimage(intensity)

    windowmaker = lambda x: np.kaiser(x, MASK_SOFTNESS)
    circles = sorted(get_circles(intensity, 3))
    virtual_order, real_order, zero_order = circles

    for circle in circles:
        print("Circle %s, %s, %s:" % circle)
        centered = get_centered(intensity, circle[1])
        showimage(centered)

    centered = get_centered(array, real_order[1])
    showimage(equalize(centered))

    window = get_holed_window(windowmaker, real_order[2] * MASK_R_SCALE, 10)
    mask = get_mask(shape, window, real_order[1])

#    window = windowmaker(zero_order[2] * MASK_R_SCALE)
#    mask *= 1 - get_mask(shape, window, zero_order[1])

    masked = get_centered(mask * array)

#    showimage(logscale(mask * 255))
#    showimage(equalize(masked))

    return masked


@Cache("pea.ifft.pickle")
def ifft(array):
    return np.fft.ifft2(array)


@Cache("pea.get_pea")
def get_pea(hologram, distance, alpha=90, beta=90):
    """
    1. hologram x ref_beam
    2. shifted_fft(1)
    3. automask(2)
    4. center(3)
    5. propagation_factor_array(M) x 5
    6. shifted_ifft(5)
    """

#    showimage(logscale(hologram))
    shape = hologram.shape
    row, col = np.ogrid[-256:256:1., -256:256:1.]
    reference = exp(1j * K * (cos(alpha) * col * DX + cos(beta) * row * DY))
    rhologram = reference * hologram
#    showimage(logscale(reference))
    frh = get_shiftedfft(rhologram)
    cfrh = get_centered(frh)
    masked = apply_mask(frh)
    phase_correction_factor = sqrt(K * 1 - (LAMBDA * 232.7920143 * row) -
        (LAMBDA * 230.8658393 * col))
    propagation_array = exp(1j * phase_correction_factor * distance)
    propagation_array = get_centered(propagation_array)
    propagated = propagation_array * masked
    reconstructed = ifft(propagated)
    wrapped = angle(reconstructed)
    return wrapped
    unwrapped = unwrap(wrapped)


def frange(start, stop, step=None, amount=None):
    assert not step or not amount

    if amount == None:
        if step == None:
            step = 1.
        amount = int((stop - start) / step)
    elif amount == 1:
        start += (stop - start) / 2.
        step = 0.
    else:
        step = (stop - start) / (amount - 1.)

    for stepn in xrange(amount):
        yield start + stepn * step


def main():
    images = [misc.imread(filename, True) for filename in sys.argv[1:]]
    for image in images:
        for alpha in frange(89.39-1.25, 89.39+1.25, amount=01):
            for beta in frange(110.57-02.5, 110.57+02.5, amount=10):
                for distance in frange(.09, .13, amount=1):
                    print("Distance: %0.2f, Alpha: %0.2fº, Beta: %0.2fº:" %
                        (distance, alpha, beta))
                    alphar = alpha / (180 / tau)
                    betar = beta / (180 / tau)
                    pea = get_pea(image, distance, alphar, betar)
                    showimage(equalize(pea))


if __name__ == "__main__":
    exit(main())
