#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from itertools import groupby, izip, count
from scipy import ndimage
import Image as pil
import ImageOps
import numpy as np
import operator


def get_centered(array, center=None, mode='wrap'):
    """
    Shift the given array to make the given point be the new center.
    If center is None the center of mass is used.
    """
    center = center or ndimage.center_of_mass(array)
    rows, cols = array.shape
    rowcc = rows / 2.
    colcc = cols / 2
    rowc, colc = center
    drows = -(rowc - rowcc)
    dcols = -(colc - colcc)
    
    centered = ndimage.shift(array, (drows, dcols), mode=mode)

#    def out2in((outrow, outcol)):
#        return outrow + drows, outcol + dcols

#    centered = ndimage.geometric_transform(array, out2in, array.shape, 
#        mode="wrap")

    return centered


def get_intensity(array):
    return array.real ** 2 + array.imag ** 2


def logscale(image):
    array = np.asarray(image, dtype=float)
    array -= array.min()
    array *= np.expm1(1) / array.max()
    array = np.log1p(array)
    array *= 255
    image = pil.fromarray(array.astype('uint8'))
    return image


#def toLmode(image):
#    image = image.convert("F")
#    array = np.asarray(image, dtype=float)
#    array -= array.min()
#    array *= 255 / array.max()
#    image = pil.fromarray(array.astype('uint8'))
#    return image


def equalizearray(array):
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
    if value:
        array *= 255. / value
    array = array.reshape(shape)
    return array


def equalize(image):
    if isinstance(image, pil.Image):
        if image.mode in ("F"):
            return equalizearray(np.asarray(image))
        elif image.mode in ("RBGA"):
            image = image.convert("RBG")
        return ImageOps.equalize(image)
    else:
        return equalizearray(image)


#def autocontrast(image):
#    if image.mode in ("F"):
#        image = toLmode(image)
#    elif image.mode in ("RBGA"):
#        image = image.convert("RBG")
#    return ImageOps.autocontrast(image)
