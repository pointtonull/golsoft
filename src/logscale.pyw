#!/usr/bin/env python
#-*- coding: UTF-8 -*-

import autopipe
import Image
import numpy as np
import sys
import ImageOps

def logscale(image):
    image = image.convert("F")
    array = np.asarray(image, dtype=float)
    array -= array.min()
    array *= np.expm1(1) / array.max()
    array = np.log1p(array)
    array *= 255
    image = Image.fromarray(array.astype('uint8'))
    return image

def fftlogscale(image):
    image = image.convert("F")
    array = np.asarray(image, dtype=float)
    array = np.fft.fft2(array)
    array -= array.min()
    array *= np.expm1(1) / array.max()
    array = np.log1p(array)
    array = np.fft.ifft2(array)
    array *= 255
    image = Image.fromarray(array.astype('uint8'))
    return image


def toLmode(image):
    image = image.convert("F")
    array = np.asarray(image, dtype=float)
    array -= array.min()
    array *= 255 / array.max()
    image = Image.fromarray(array.astype('uint8'))
    return image


def equalizefloat(image):
    image = image.convert("F")
    array = np.asarray(image)
    shape = array.shape
    array = array.flatten()
    sorters = array.argsort()
    for order, place in enumerate(sorters):
        array[place] = order
    array *= 255. / order
    print(array.min(), array.max())
    array = array.reshape(shape)
    image = Image.fromarray(array.astype("uint8"))
    return image


def main():
    image = Image.open(sys.argv[1])

    print("Original image:")
    autopipe.showimage(image)

    print("On gray tones:")
    autopipe.showimage(image.convert("F"))

    print("On log-scaled gray tones:")
    autopipe.showimage(logscale(image))

    print("With auto-contrast:")
    autopipe.showimage(ImageOps.autocontrast(toLmode(image)))

    print("With equalized histogram:")
    autopipe.showimage(ImageOps.equalize(toLmode(image)))

    print("With float equalized histogram:")
    autopipe.showimage(equalizefloat(image))

if __name__ == "__main__":
    exit(main())
