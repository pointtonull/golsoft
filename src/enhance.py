#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from itertools import groupby, izip, count
import Image as pil
import ImageOps
import numpy as np
import operator


def logscale(image):
    array = np.asarray(image, dtype=float)
    array -= array.min()
    array *= np.expm1(1) / array.max()
    array = np.log1p(array)
    array *= 255
    image = pil.fromarray(array.astype('uint8'))
    return image


def toLmode(image):
    image = image.convert("F")
    array = np.asarray(image, dtype=float)
    array -= array.min()
    array *= 255 / array.max()
    image = pil.fromarray(array.astype('uint8'))
    return image


def equalizefloat(array):
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
#    print "frequencies: %d" % value
    return array


def equalize(image):
    if isinstance(image, pil.Image):
        if image.mode in ("F"):
            return equalizefloat(np.asarray(image))
        elif image.mode in ("RBGA"):
            image = image.convert("RBG")
        return ImageOps.equalize(image)
    else:
        return equalizefloat(image)

def autocontrast(image):
    if image.mode in ("F"):
        image = toLmode(image)
    elif image.mode in ("RBGA"):
        image = image.convert("RBG")
    return ImageOps.autocontrast(image)

if __name__ == "__main__":
    from mayavi.api import Engine
    from mayavi.sources.api import ArraySource
    from mayavi.filters.api import WarpScalar, PolyDataNormals
    from mayavi.modules.api import Surface
    import scipy
    array = logscale(scipy.lena())
    print max(array)
    engine = Engine()
    engine.start()
    scene = engine.new_scene()
    source = ArraySource(scalar_data=array)
    engine.add_source(source)
    warp = WarpScalar()
    engine.add_filter(warp, obj=source)
    normals = PolyDataNormals()
    engine.add_filter(normals, obj=warp)
    surf = Surface()
    engine.add_module(surf, obj=normals)
