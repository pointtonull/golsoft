#!/usr/bin/env python
#-*- coding: UTF-8 -*-

import autopipe
import Image
import numpy
import sys


def logscale(image):
    array = numpy.asarray(image.convert("F"), dtype=float)
    array -= array.min()
    array *= (numpy.exp(1) - 1) / array.max()
    array = numpy.log1p(array)
    array *= 255
    image = Image.fromarray(array.astype('uint8'))
    return image

def main():
    image = Image.open(sys.argv[1])
    print("Original image:")
    autopipe.showimage(image)
    print("On gray tones:")
    autopipe.showimage(image.convert("F"))
    image = logscale(image)
    print("On log-scaled gray tones:")
    autopipe.showimage(image)

if __name__ == "__main__":
    exit(main())
