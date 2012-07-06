#!/usr/bin/env python
#-*- coding: UTF-8 -*-

import sys

import numpy as np
import matplotlib.pyplot as plt

from autopipe import showimage
from dft import get_shifted_dft
from image import imread
from pea import get_auto_mask, get_propagation_array

def get_phase_diff2(spectrum):
    # input of number of range bins nrgb
    pass

    # generate sidelobe weights subwts
    pass

    # input normalization choice
    pass

    # call phase difference suboutine
    pass

    # define sub-array length nsub
    # define subarray center LI
    # unpack array x from first half of range bin
    # unpack array y grom second half of range bin
    # unpack array z from middle half of range bin
    x, y, z = partition(spectrum)

    # compute conjugate of array x
    # compute conjugate o array z
    x_conjugate = x.conjugate()
    z_conjugate = z.conjugate()

    # sum A1 = array y x array x*
    # sum A2 = array z x array x*
    # sum A3 = array y x array z*
    # apply weights
    # alternate sign
    sum1 = -sum(y * x_conjugate, 0)
    sum2 = -sum(z * x_conjugate, 0)
    sum3 = -sum(y * z_conjugate, 0)

    # sum magnitude of array in sum
    magsum1 = np.abs(sum1)
    magsum2 = np.abs(sum2)
    magsum3 = np.abs(sum3)
    # zero fill to fft size
    # call fft subroutine
    # calculate fft of suma1 suma2 suma3
    dftsum1 = get_shifted_dft(sum1)
    dftsum2 = get_shifted_dft(sum2)
    dftsum3 = get_shifted_dft(sum3)
    # magnitude of filter output
    magdft1 = np.abs(dftsum1)
    magdft2 = np.abs(dftsum2)
    magdft3 = np.abs(dftsum3)
    # if normalize:
    #     divide por insum
    mag1 = magsum1 / sum(magsum1)
    mag2 = magsum2 / sum(magsum1)
    mag3 = magsum3 / sum(magsum1)

#    #experimental
#    mag2[mag2 > (mag2.min() + mag2.ptp() / 2.)] = 0
#    mag3[mag3 > (mag3.min() + mag3.ptp() / 2.)] = 0

    convolution = np.convolve(mag2, mag3)

    figure = plt.figure()

#    plt.cla()
#    plt.plot(mag1, c="blue")
#    showimage(figure)

#    plt.cla()
#    plt.plot(mag2, c="green")
#    showimage(figure)

#    plt.cla()
#    plt.plot(mag3, c="red")
#    showimage(figure)

    plt.cla()
    plt.plot(mag2 / mag3, c="red")
    showimage(figure)

#    plt.cla()
#    plt.plot(convolution)
#    showimage(figure)

    return convolution.argmax()
    # sum over range bin
    # call peak detection routine
    # locate maximun index
    # interpolate
    # calculate peak index
    # calculate tau
    # calculate cuadratic error
    # calculate cubic error


def partition(spectrum):
    center = np.array(spectrum.shape) / 2.
    x = spectrum[:center[1], :]
    y = spectrum[center[1]:, :]
    z = spectrum[center[1] / 2:center[1] / 2 + center[1], :]
    return x, y, z


def peak_detect(left, rigth):
    s32 = x.transpose()
    left = left.transpose()

    s34 = s32 * z
    product = left * rigth

    s36 = adaptative_scaling(s34)
    scaled_product = adaptative_scaling(product)

    s38 = wts(s36)
    print s38

    s40 = - s38
    s42 = zero_fill(s40)
    print s42
    s44 = get_shifted_dft(s42)
    s46 = mag(s44)
    s48 = mag(s40)
    s50 = 1/sum(s48)
    s52 = x(s46, s50)
    s54 = sum(s52)
    s56 = peak_detect(s54)


def get_phase_diff(spectrum):

    azimuts = spectrum.sum(1)
    x, y, z = partition(spectrum)

    def adaptative_scaling(arg):
        return arg

    def wts(arg):
        return arg

    def zero_fill(arg):
        return arg

    s58 = z.transpose()
#    s60 = (...)
    s62 = peak_detect(x, y)

    s66 = quad_cubic_error_estimates(s56, s62)


def main():
    from image import equalize
    images = [(filename, imread(filename, True)) for filename in sys.argv[1:]]

    for filename, image in images:
        shape = image.shape
        print("\nOriginal image: %s" % filename)
        spectrum = get_shifted_dft(image)
        mask, masked_spectrum, centered = get_auto_mask(spectrum,
            softness=0, radious_scale=1.5)
        showimage(image, equalize(masked_spectrum))
        for dist in range(-20, 20 + 1):
            dist = dist / 2.
            propagated = masked_spectrum * get_propagation_array(shape, dist)
            phase_diff = get_phase_diff2(propagated)
            print("%6.2f  %5.2f" % (dist, phase_diff))

if __name__ == "__main__":
    exit(main())
