#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from itertools import groupby, izip, count
from scipy import ndimage, misc
import Image as pil
import ImageOps
import cache
import numpy as np
import operator

VERBOSE = 0


def open_raw(filename, aspectratio=1):
    bits = open("../docs/holos/32h3_400.raw").read()
    length = len(bits)
    cols = int(round((length * aspectratio) ** .5))
    rows = length / cols
    if length != cols * rows:
        raise ValueError("incorrect aspectratio")
    array = np.array([ord(char) for char in bits])
    array = array.reshape((rows, cols))
    return array


def imread(filename, flatten=True, aspectratio=1):
    try:
        array = misc.imread(filename, flatten)
    except IOError:
        array = open_raw(filename, aspectratio)
    return array


@cache.hybrid
def get_centered(array, center=None, mode='wrap'):
    """
    Shift the given array to make the given point be the new center.
    If center is None the center of mass is used.
    mode can be 'constant', 'nearest', 'reflect' or 'wrap'.
    """

    if center:
        rows, cols = array.shape
        rowcc = rows / 2.
        colcc = cols / 2
        rowc, colc = center
        drows = rowcc - rowc
        dcols = colcc - colc
        shift = (drows, dcols)

    else:
        if issubclass(array.dtype.type, complex):
            intensity = get_intensity(array)
            shift = get_shift_to_center_of_mass(intensity, mode)
        else:
            shift = get_shift_to_center_of_mass(array, mode)

    if issubclass(array.dtype.type, complex):
        real = ndimage.shift(array.real, shift, mode=mode)
        imag = ndimage.shift(array.imag, shift, mode=mode)
        centered = real + 1j * imag
    else:
        centered = ndimage.shift(array, shift, mode=mode)

    return centered


@cache.hybrid
def get_shift_to_center_of_mass(array, mode="wrap"):
    """
    Calcules the shift of the center of mass relative to the center of the image
    """
    if array.ndim > 1:
        shift = [get_shift_to_center_of_mass(array.sum(dim))
            for dim in range(array.ndim)][::-1]
        return shift
    else:
        center = array.shape[0] / 2.
        total_shift = 0
        centered = array
        for step in xrange(100):
            center_of_mass = ndimage.center_of_mass(centered)
            shift = center - center_of_mass[0]
            eshift = shift * 2 ** .5
            if abs(eshift) < 1:
                break
            total_shift += eshift
            centered = ndimage.shift(centered, eshift, mode=mode)
        return total_shift


def get_intensity(array):
    return array.real ** 2 + array.imag ** 2


def logscale(array):
    if issubclass(array.dtype.type, complex):
        array = get_intensity(array)
    array = array.astype(float)
    array -= array.min()
    array *= np.expm1(1) / array.max()
    array = np.log1p(array)
    array *= 255.
    return array


def equalizearray(array):
    if issubclass(array.dtype.type, complex):
        array = get_intensity(array)
    array = array.astype(float)
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


@cache.hybrid
def equalize(image):
    if isinstance(image, pil.Image):
        if image.mode in ("F"):
            return equalizearray(np.asarray(image))
        elif image.mode in ("RBGA"):
            image = image.convert("RBG")
        return ImageOps.equalize(image)
    else:
        return equalizearray(image)
