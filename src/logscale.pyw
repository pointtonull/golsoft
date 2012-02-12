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
    print(array.min(), array.max())

    array = np.fft.fft2(array)
    print(array.min(), array.max())

    array -= array.min()
    print(array.min(), array.max())

    array *= np.expm1(1) / array.max()
    print(array.min(), array.max())

    array = np.log1p(array)
    print(array.min(), array.max())

    array = np.fft.ifft2(array)
    print(array.min(), array.max())

    array *= 255
    print(array.min(), array.max())

    image = Image.fromarray(array.astype('uint8'))
    return image

def toLmode(image):
    image = image.convert("F")
    array = np.asarray(image, dtype=float)
    array -= array.min()
#    array *= 255 / array.max()
    image = Image.fromarray(array.astype('uint8'))
    return image


def main():
    image = Image.open(sys.argv[1])

#    print("Original image:")
#    autopipe.showimage(image)

    print("On gray tones:")
    autopipe.showimage(image.convert("F"))

#    print("On log-scaled gray tones:")
#    autopipe.showimage(logscale(image))

    print("With auto-contrast:")
    autopipe.showimage(ImageOps.autocontrast(toLmode(image)))

    print("With equalized histogram:")
    autopipe.showimage(ImageOps.equalize(toLmode(image)))

if __name__ == "__main__":
    exit(main())
