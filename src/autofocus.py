#!/usr/bin/env python
#-*- coding: UTF-8 -*-


def partition(spectrum):
    x = spectrum[:center, :]
    y = spectrum[center:, :]
    z = spectrum[center/2:center/2+center, :]
    return x, y, z


def get_phase_diff(spectrum):

    # input of number of range bins nrgb
    # generate sidelobe weights subwts
    # input normalization choice
    # call phase difference suboutine
    # define sub-array length nsub
    # define subarray center LI
    # unpack array x from first half of range bin
    # unpack array y grom second half of range bin
    # unpack array z from middle half of range bin
    # compute conjugate of array x
    # compute conjugate o array z
    # sum A1 = array y x array x*
    # sum A2 = array Z x array x*
    # sum A3 = array y x array z*
    # apply weights
    # alternate sign
    # sum magnitude of array in sum
    # zero fill to fft size
    # call fft subroutine
    # calculate fft of suma1 suma2 suma3
    # magnitude of filter output
    # if normalize:
    #     divide por insum
    # sum over range bin
    # call peak detection routine
    # locate maximun index
    # interpolate
    # calculate peak index
    # calculate tau
    # calculate cuadratic error
    # calculate cubic error

    azimuts = spectro.sum(1)
    center = azimuts.shape / 2

    x = azimuts[:center]
    y = azimuts[center:]
    z = azimuts[center - center/2:center + center/2]

    def peak_detect(left, rigth):
        s32 = x.transpose()
        left = left.transpose()

        s34 = s32 * z
        product = left * rigth

        s36 = adaptative_scaling(34)
        scaled_product = adaptative_scaling(product)

        s38 = wts(36)

        s40 = +-+(38)
        s42 = zero_fill(40)
        s44 = fft(42)
        s46 = mag(44)
        s48 = mag(40)
        s50 = 1/sum(48)
        s52 = x(46, 50)
        s54 = sum(52)
        s56 = peak_detect(54)

    s58 = z.transpose()
    s60 = (...)
    s62 = peak_detect()

    s66 = quad_cubic_error_estimates(56, 62)


def main():
    image = None

    spectro = fft(image)
    phase_diff = get_phase_diff(spectro)
    azzimut = fft(spectro (*) phase_diff)


if __name__ == "__main__":
    exit(main())
