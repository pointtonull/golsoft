#!/usr/bin/env python
#-*- coding: UTF-8 -*-

"""
Simple implementation of the Angular Spectrum Method to reconstruct lensless
holograms
"""

from numpy import exp
import numpy as np

from autofocus import guess_focus_distance
from automask import get_circles, get_auto_mask
from dependences import Datum, Depends
from dft import get_shifted_dft, get_idft, get_shifted_idft
from propagation import get_propagation_array
from image import get_intensity, imread
from unwrap import unwrap_wls
import cache



tau = 6.283185307179586
LAMBDA = 6.328e-07 # wave length
DX = 8.39e-6
DY = 8.46e-6
K = tau / LAMBDA # wave number
EPSILON = 1e-16


def angle2(array):
    raise DeprecationWarning("use get_phase insteat")
    return np.arctan2(array.real, array.imag)


def get_phase(array):
    return np.arctan2(array.real, array.imag)


def get_module(array):
    return np.abs(array)




def get_refbeam(shape, cos_alpha=EPSILON, cos_beta=EPSILON):
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
    ref_beam = get_refbeam(shape, cos_alpha, cos_beta)
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
 

def get_distance(point1, point2):
    distance = ((point1[0] - point2[0]) ** 2
        + (point1[1] - point2[1]) ** 2) ** .5
    return distance


def get_peak_coords(spectrum):
    # TODO: must select by position not by intensity
    shape = spectrum.shape
    center = [dim / 2. for dim in shape]

    intensity = get_intensity(spectrum)
    circles = [(-get_distance(center, circle[1]), circle[0], circle[1])
        for circle in get_circles(intensity, 2, 20)]
    circles.sort()
    peak = circles[0][2]
    peaks_row = (peak[0] - center[0]) / float(center[0])
    peaks_col = (peak[1] - center[1]) / float(center[1])
    return peaks_row, peaks_col


@cache.hybrid(reset=0)
def calculate_director_cosines(spectrum):
    """
    Calculate the director cosines using the spectral proyection formula
    """
    peak = get_peak_coords(spectrum)
    freq_rows, freq_cols = peak
    freq_rows /= 2 * DY
    freq_cols /= 2 * DX
    cos_alpha = freq_cols * LAMBDA
    cos_beta = freq_rows * LAMBDA

    return cos_alpha, cos_beta



class PEA(object):

    def __init__(self, filename=None):
        if filename:
            self.filename = filename

    filename = Datum()
    @Depends(filename)
    def image(self):
        print("Loading image")
        return imread(self.filename, True)


    @Depends(image)
    def ispectrum(self):
        print("DFT(image)")
        return get_shifted_dft(self.image)


    use_autocosines = Datum(True)
    user_cosines = Datum((0, 0))
    @Depends(ispectrum, use_autocosines, user_cosines)
    def cosines(self):
        print("Director cosines")
        if self.use_autocosines:
            return calculate_director_cosines(self.ispectrum)
        else:
            return self.user_cosines


    @Depends(image, cosines)
    def refbeam(self):
        print("Calculating refbeam")
        return get_refbeam(self.image.shape, *self.cosines)


    use_refbeam = Datum(False)
    @Depends(image, refbeam, use_refbeam)
    def hologram(self):
        print("R-Hologram")
        if self.use_refbeam:
            return self.refbeam * self.image
        else:
            return self.image


    @Depends(image, cosines)
    def spectrum(self):
        print("DFT(R-Hologram)")
        if self.use_refbeam:
            return get_shifted_dft(self.hologram)
        else:
            return self.ispectrum


    softness = Datum(0)
    order_scale = Datum(0.8)
    use_zeromask = Datum(True)
    zero_scale = Datum(1.2)
    use_cuttop = Datum(False)
    cuttop = Datum(0.005)
    @Depends(spectrum, order_scale, use_zeromask, zero_scale, softness,
        use_cuttop, cuttop)
    def masking(self):
        print("Masking")
        zero_scale = self.zero_scale if self.use_zeromask else 0
        cuttop = self.cuttop if self.use_cuttop else 0
        mask, masked, centered = get_auto_mask(self.spectrum,
            self.softness, self.order_scale, zero_scale, cuttop)
        return mask, masked, centered


    @Depends(masking)
    def mask(self):
        print("Mask")
        return self.masking[0]


    @Depends(masking)
    def masked_spectrum(self):
        print("Masked spectrum")
        return self.masking[1]


    @Depends(masking)
    def centered_spectrum(self):
        print("Centered spectrum")
        return self.masking[2]


    @Depends(ispectrum)
    def auto_distance(self):
        print("Auto distance")
        mask, masked, centered = get_auto_mask(self.ispectrum)
        distance = guess_focus_distance(masked)
        return distance


    use_autofocus = Datum(True)
    user_distance = Datum(0.05)
    @Depends(use_autofocus, user_distance, auto_distance)
    def distance(self):
        print("Distance")
        if self.use_autofocus:
            return self.auto_distance
        else:
            return self.user_distance


    @Depends(spectrum, distance)
    def propagation(self):
        print("Propagation")
        return get_propagation_array(self.spectrum.shape, self.distance)


    @Depends(masked_spectrum, propagation)
    def propagated(self):
        print("Propagated")
        return self.masked_spectrum * self.propagation


    @Depends(propagated)
    def reconstructed(self):
        print("IDFT(Propagated)")
        return get_shifted_idft(self.propagated)

    @Depends(reconstructed)
    def module(self):
        print("Module")
        return get_module(self.reconstructed)

    @Depends(reconstructed)
    def phase(self):
        print("Phase")
        return get_phase(self.reconstructed)

    unwrapper = Datum(unwrap_wls)
    @Depends(phase, module, unwrapper)
    def unwrapped_phase(self):
        print("Unwrapped phase")
        return self.unwrapper(self.phase, self.module)
