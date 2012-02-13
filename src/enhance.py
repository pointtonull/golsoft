#!/usr/bin/env python
#-*- coding: UTF-8 -*-

import Image
import ImageOps
from itertools import groupby, izip, count
import operator
import numpy as np


def logscale(image):
    array = np.asarray(image, dtype=float)
    array -= array.min()
    array *= np.expm1(1) / array.max()
    array = np.log1p(array)
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
    array = np.asarray(image)
    shape = array.shape
    array = array.flatten()
    sorters = array.argsort()
    array.sort()
    zippeds = izip(array, sorters)
    groups = groupby(zippeds, operator.itemgetter(0))
    counter = count()
    for ovalue, group in groups:
        value = counter.next()
        for ovalue, pos in list(group):
            array[pos] = value
    array *= 255. / value
    array = array.reshape(shape)
    image = Image.fromarray(array.astype("uint8"))
    print "frequencies: %d" % value
    return image


def equalize(image):
    if image.mode in ("F"):
        return equalizefloat(image)
    elif image.mode in ("RBGA"):
        image = image.convert("RBG")
    return ImageOps.equalize(image)


def autocontrast(image):
    if image.mode in ("F"):
        image = toLmode(image)
    elif image.mode in ("RBGA"):
        image = image.convert("RBG")
    return ImageOps.autocontrast(image)

