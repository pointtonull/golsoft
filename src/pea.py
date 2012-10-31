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
from dft import get_shifted_dft, get_shifted_idft, get_module, get_phase
from image import get_intensity, imread, subtract, limit_size, equalize
from image import phase_denoise
from propagation import get_propagation_array
from unwrap import unwrap_wls
import cache


tau = 6.283185307179586 # two times sexier than pi


def get_refbeam(shape, cos_alpha, cos_beta, wavelength, (dx, dy)):
    """
    Generate a reference beam array given the shape of the hologram and the
    directors angles
    """
    wavenumber = tau / wavelength
    maxrow = shape[0] / 2
    maxcol = shape[1] / 2
    minrow, mincol = -maxrow, -maxcol
    row, col = np.ogrid[minrow:maxrow:1., mincol:maxcol:1.]

    ref_beam = exp(1j * wavenumber * (cos_alpha * col * dx + cos_beta *
        row * dy))
    
    return ref_beam


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
def calculate_director_cosines(spectrum, wavelength, (dx, dy)):
    """
    Calculate the director cosines using the spectral proyection formula
    """
    peak = get_peak_coords(spectrum)
    freq_rows, freq_cols = peak
    freq_rows /= 2 * dy
    freq_cols /= 2 * dx
    cos_alpha = freq_cols * wavelength
    cos_beta = freq_rows * wavelength

    return cos_alpha, cos_beta



class PEA(object):

    def __init__(self, filename=None):

        if filename:
            self.filename_holo = filename


    resolution_limit = Datum(1)
    filename_holo = Datum()
    @Depends(filename_holo, resolution_limit)
    def image_holo(self):
        print("Loading hologram image")
        image = imread(self.filename_holo, True)
        image = limit_size(image, self.resolution_limit)
        return image


    filename_ref = Datum("")
    @Depends(filename_ref, resolution_limit)
    def image_ref(self):
        print("Loading reference image")
        image = imread(self.filename_ref, True)
        image = limit_size(image, self.resolution_limit)
        return image


    filename_obj = Datum("")
    @Depends(filename_obj, resolution_limit)
    def image_obj(self):
        print("Loading object image")
        image = imread(self.filename_obj, True)
        image = limit_size(image, self.resolution_limit)
        return image


    equalize_image = Datum(True)
    @Depends(image_holo, image_ref, image_obj, equalize_image)
    def image(self):
        print("Calculating hologram")
        image = self.image_holo
        if self.filename_ref:
            image = subtract(image, self.image_ref)
        if self.filename_obj:
            image = subtract(image, self.image_obj)
        if self.equalize_image:
            image = equalize(image)
        return image


    @Depends(image)
    def ispectrum(self):
        print("DFT(image)")
        return get_shifted_dft(self.image)


    use_autocosines = Datum(True)
    user_cosines = Datum((0, 0))
    wavelength = Datum(650e-9)
    dx = Datum(3e-6)
    dy = Datum(3e-6)
    @Depends(ispectrum, use_autocosines, user_cosines, wavelength, dx, dy)
    def cosines(self):
        print("Director cosines")
        if self.use_autocosines:
            return calculate_director_cosines(self.ispectrum,
                self.wavelength, (self.dx, self.dy))
        else:
            return self.user_cosines


    @Depends(image, cosines, wavelength)
    def refbeam(self):
        print("Calculating refbeam")
        cos_alpha, cos_beta = self.cosines
        refbeam = get_refbeam(self.image.shape, cos_alpha, cos_beta,
            self.wavelength, (self.dx, self.dy))
        return refbeam


    use_refbeam = Datum(False)
    @Depends(image, refbeam, use_refbeam)
    def hologram(self):
        print("R-Hologram")
        if self.use_refbeam:
            return self.refbeam * self.image
        else:
            return self.image


    @Depends(hologram, cosines)
    def spectrum(self):
        print("Spectrum (dft(hologram))")
        if self.use_refbeam:
            return get_shifted_dft(self.hologram)
        else:
            return self.ispectrum


    mask_softness = Datum(0)
    mask_order_scale = Datum(0.8)
    mask_use_zeromask = Datum(True)
    mask_zero_scale = Datum(1.2)
    mask_use_cuttop = Datum(False)
    mask_cuttop = Datum(0.005)
    @Depends(spectrum, mask_order_scale, mask_use_zeromask, mask_zero_scale,
        mask_softness, mask_use_cuttop, mask_cuttop)
    def masking(self):
        print("Masking")
        zero_scale = self.mask_zero_scale if self.mask_use_zeromask else 0
        cuttop = self.mask_cuttop if self.mask_use_cuttop else 0
        mask, masked, centered = get_auto_mask(self.spectrum,
            self.mask_softness, self.mask_order_scale, zero_scale, cuttop)
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


    @Depends(ispectrum, wavelength, dx, dy)
    def auto_distance(self):
        print("Auto distance")
        mask, masked, centered = get_auto_mask(self.ispectrum)
        distance = guess_focus_distance(masked, self.wavelength,
            (self.dx, self.dy))
        return distance


    use_autofocus = Datum(False)
    user_distance = Datum(0)
    @Depends(use_autofocus, user_distance, auto_distance)
    def distance(self):
        print("Distance")
        if self.use_autofocus:
            return self.auto_distance
        else:
            return self.user_distance


    @Depends(spectrum, distance, wavelength, dx, dy)
    def propagation(self):
        print("Propagation")
        return get_propagation_array(self.spectrum.shape, self.distance,
            self.wavelength, (self.dx, self.dy))


    propagate = Datum(True)
    @Depends(propagate, masked_spectrum, propagation)
    def propagated(self):
        if self.propagate:
            print("Propagated")
            return self.masked_spectrum * self.propagation
        else:
            return self.masked_spectrum


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


    phase_denoise = Datum(0)
    @Depends(phase, phase_denoise)
    def denoised_phase(self):
        if self.phase_denoise:
            print("Denoising phase")
            return phase_denoise(self.phase, self.phase_denoise)
        else:
            return self.phase


    unwrapper = Datum(unwrap_wls)
    @Depends(phase, module, unwrapper, phase_denoise)
    def unwrapped_phase(self):
        print("Unwrapping phase")
        if self.unwrapper.func_code.co_argcount == 1:
            return self.unwrapper(self.denoised_phase)
        else:
            return self.unwrapper(self.denoised_phase, self.module)
