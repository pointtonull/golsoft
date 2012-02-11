#!/usr/bin/env python
#-*- coding: UTF-8 -*-

#import autopipe
import Image
import numpy
import sys


def logscale(image):
    array = numpy.array(image, dtype=float)
    print array.min(), array.max()
    array -= array.min()
    print array.min(), array.max()
    array *= numpy.exp(1) / array.max()
    print array.min(), array.max()
    array = numpy.log1p(array)
    print array.min(), array.max()
    array *= 255
    print array.min(), array.max()
    image = Image.fromarray(array)
    return image

def main():
    image = Image.open(sys.argv[1])
    image.show("Original")
    image = logscale(image)
    image.show("Logscaled")
#    autopipe.pipe.writeimage(image)

if __name__ == "__main__":
    exit(main())
