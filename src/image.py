#!/usr/bin/env python
#-*- coding: UTF-8 -*-

import operator

from itertools import groupby, izip, count

import Image as pil
import ImageOps
import gdal
import numpy as np

from numpy import sin, cos, exp, log, arctan2
from scipy import misc, ndimage
from scipy.misc import imresize
from scipy.ndimage import geometric_transform
from skimage.filter import canny

from minimize import generic_minimizer
import cache

VERBOSE = 0
tau = np.pi * 2 # two times sexier than pi

#TODO:
#  * implement:
#    - imshow   as in autopipe
#    - imsave


def correct_cmos_stripes(phase, position=0.5, denoise=5):
    """
    Al capturar hologramas con CMOS pueden aparecer franjas horizontales que
    deforman la phase. Estas franjas pueden ser aisladas por su regularidad
    horizotal.

    position: medida de posición no central, puede ser cualquier flotante entre 0 y 1:
        1.0 es el máximo
        0.5 (default) equivale la mediana, al segundo cuantil y al 50º percentil
        0.0 es el mínimo

    denoise: tamaño del filtro gausiano aplicado a los valores de las franjas
    """
    assert 0 <= position <= 1
    values = np.sort(phase)
    denoised = ndimage.filters.gaussian_filter(values, denoise)
    rows, cols = phase.shape
    col = cols * position
    strip = denoised[:, col].reshape((rows, 1))
    return phase - strip


def choice_density_map(density_map, random=None):
    cumsum = np.cumsum(density_map)
    if random is None:
        random = np.random.rand()
    pos = np.searchsorted(cumsum, random * cumsum[-1])
    row, col = np.unravel_index(pos, density_map.shape)
    print("(%d, %d)" % (row, col))
    return row, col


def phase_denoise(phase, size=1, filter_func=ndimage.filters.median_filter):
    """
    Cuadratic denoise. Is a median filter applied on the angular space.
    """
    if size == 0:
        return phase % tau
    else:
        y_over = sin(phase)
        x_over = cos(phase)
        y_over = filter_func(y_over, size)
        x_over = filter_func(x_over, size)
        denoised = arctan2(y_over, x_over)
        
    return denoised


def phase_filter(phase, filter_func=np.var, *extra_args):
    y_over = (sin(phase) + 1) / 2.
    print y_over.ptp(), y_over.min()
    y_over = filter_func(y_over, *extra_args)
    return y_over


def auto_canny(array, average=0.15, gaussian_sigma=1, strongness=2.5):
    array -= array.min()
    array /= array.max()

    def canny_average(hard_threshold):
        soft_threshold = hard_threshold / strongness
        edges = canny(array, gaussian_sigma, hard_threshold, soft_threshold)
        return edges.mean()

    hard_threshold = 0.4
    current_average = canny_average(hard_threshold)
    epsilon = 0.001
    for iteration in xrange(50):
        if abs(current_average - average) < epsilon:
            break
        if current_average < average:
            hard_threshold /= 1.2
        else:
            hard_threshold *= 1.2
        print(hard_threshold, current_average)

    soft_threshold = hard_threshold / strongness
    return canny(array, gaussian_sigma, hard_threshold, soft_threshold)



@cache.hybrid
def get_subtract_paramns(left, right):
    """
    Returns k that minimizes:

        var(left - k * right)
    """

    def diference((k)):
        return (left - k * right).var()

    best_k = float(generic_minimizer(diference, 1))
    return best_k


def subtract(left, right):
    """
    Will operate
        left - k * right + l

    Where k and l are the values that minimizes the result.
    """

    if right is None:
        return left
    else:
        best_k = get_subtract_paramns(left, right)
        print("Subtraction left - %f * right" % best_k)
        result = left - best_k * right

        return result


def limit_size(image, limit, avoidodds=True):
    """
    Image is a numpy array.
    Resulution is a quantity:
        in pixels if >= 1000
        in megapixels if < 1000
    """

    if limit < 1000:
        limit *= 1e6

    relsize = (image.size / limit) ** -.5
    if relsize <= 1:
        new_shape = [int(round(res * relsize))
            for res in image.shape]

        if avoidodds:
            new_shape = tuple([int(res + res % 2)
                for res in new_shape])

        image = imresize(image, new_shape, 'bicubic')
        image = np.float32(image)
    return image


def get_logpolar(array, interpolation=0, reverse=False):
    """
    Returns a new array with the logpolar transfamation of array.
    Interpolation can be:
        0 Near
        1 Linear
        2 Bilineal
        3 Cubic
        4
        5
    """
    assert interpolation in range(6)
    rows, cols = array.shape
    row0 = rows / 2.
    col0 = cols / 2.
    theta_scalar = tau / cols
    max_radius = (row0 ** 2 + col0 ** 2) ** .5
    rho_scalar = log(max_radius) / cols

    def cart2logpol(dst_coords):
        theta, rho = dst_coords
        rho = exp(rho * rho_scalar)
        theta = np.pi / 2 - theta * theta_scalar
        row_from = rho * cos(theta) + row0
        col_from = rho * sin(theta) + col0
        return row_from, col_from

    def logpol2cart(dst_coords):
        xindex, yindex = dst_coords
        x = xindex - col0
        y = yindex - row0

        r = np.log(np.sqrt(x ** 2 + y ** 2)) / rho_scalar
        theta = np.arctan2(y, x)
        theta_index = np.round((theta + np.pi) * cols / tau)
        return theta_index, r

    trans = logpol2cart if reverse else cart2logpol

    logpolar = geometric_transform(array, trans, array.shape,
        order=interpolation)

    return logpolar


def get_polar(array, interpolation=0, reverse=False):
    """
    Returns a new array with the logpolar transfamation of array.
    Interpolation can be:
        0 Near
        1 Linear
        2 Bilineal
        3 Cubic
        4
        5
    """
    assert interpolation in range(6)
    rows, cols = array.shape
    row0 = rows / 2.
    col0 = cols / 2.
    theta_scalar = tau / cols
    max_radius = (row0 ** 2 + col0 ** 2) ** .5
    rho_scalar = max_radius / cols

    def cart2pol(dst_coords):
        theta, rho = dst_coords
        rho = rho * rho_scalar
        theta = np.pi / 2 - theta * theta_scalar
        row_from = rho * cos(theta) + row0
        col_from = rho * sin(theta) + col0
        return row_from, col_from

    def pol2cart(dst_coords):
        xindex, yindex = dst_coords
        x = xindex - col0
        y = yindex - row0

        r = np.sqrt(x ** 2 + y ** 2) / rho_scalar
        theta = np.arctan2(y, x)
        theta_index = np.round((theta + np.pi) * cols / tau)
        return theta_index, r

    trans = pol2cart if reverse else cart2pol

    polar = geometric_transform(array, trans, array.shape,
        order=interpolation)

    return polar


def open_raw(filename):
    known_resolutions = {
        5038848: (1944, 2592, "bayer"),
        262144: (512, 512, "mono"),
        266638: (512, 520, "mono"),
    }

    bits = open(filename, "rb").read()
    length = len(bits)

    if length in known_resolutions:
        rows, cols, method = known_resolutions[length]
        array = np.array([ord(char) for char in bits])
        array = array.reshape((rows, cols))

        if method == "bayer":
            #TODO: implement Malvar-He-Cutler Bayer demosaicing
            print("Identified %s as bayer raw." % filename)
            array0 = array[0::2, 0::2]
            array1 = array[0::2, 1::2]
            array2 = array[1::2, 0::2]
            array3 = array[1::2, 1::2]
            red = array1
            green = (array0 + array3) / 2
            blue = array2
            array = np.array([red, green, blue])

        return array

    else:
        raise IOError("unknown resolution on raw file %s (%d pixels)" %
            (filename, length))


def open_gdal(filename):
    dataset = gdal.Open(filename)
    array = dataset.ReadAsArray()
    return array


def imread(filename, flatten=True):
    if filename.endswith(".raw"):
        array = open_raw(filename)
    else:
        try:
            array = misc.imread(filename, flatten)
        except IOError:
            array = open_gdal(filename)
            array = array[:3, :, :] # alpha shift
            if flatten:
                array = array.mean(0)
    return array


def derotate(array):
    rows, cols = array.shape
    polar_array = get_logpolar(array, 0)
    rows_sum = polar_array.sum(1)
    maxcol = - rows_sum.argmax()
    rows_sum = ndimage.shift(rows_sum, maxcol, order=0, mode="wrap")
    rows_shift = maxcol
    angle = (-360. * rows_shift) / rows
    derotated = ndimage.rotate(array, angle, reshape=False)
    return derotated


def evenshape(array, shrink=False):
    if not shrink:
        newshape = [dim + 1 - dim % 2 for dim in array.shape]
        newarray = np.zeros(newshape)
        newarray[:array.shape[0],:array.shape[1]] = array
    else:
        newshape = [dim - 1 + dim % 2 for dim in array.shape]
        newarray = array[:newshape[0], :newshape[1]]

    return newarray



def imwrite(array, filename):
    return misc.imsave(filename, array)


def get_centered(array, center=None, mode='wrap', reverse=False):
    """
    Shift the given array to make the given point be the new center.
    If center is None the center of mass is used.
    mode can be 'constant', 'nearest', 'reflect' or 'wrap'.
    
    inverse False:  center -> current_center
    inverse True:   current_center -> center
    
    """

    if center:
        rows, cols = array.shape
        rowcc = int(round(rows / 2.))
        colcc = int(round(cols / 2.))
        rowc, colc = center
        if reverse:
            drows = rowc - rowcc
            dcols = colc - colcc
        else:
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

        shift = int(round(total_shift))

        return shift


def get_intensity(array):
    return array.real ** 2 + array.imag ** 2


def logscale(array):
    array = array.copy()
    if issubclass(array.dtype.type, complex):
        array = get_intensity(array)
    array = array.astype(float)
    array -= array.min()
    array *= np.expm1(1) / array.max()
    array = np.log1p(array)
    array *= 255.
    return array


def normalize(array):
    """
    Apply linears tranformations to ensure all the values are in [0, 255]
    """
    array = array.copy()
    if issubclass(array.dtype.type, complex):
        array = get_intensity(array)
    array -= array.min()
    array *= 255. / array.max()
    return array


def equalizearray(array):
    """
    Equalize the array histogram
    """
    array = normalize(array)
    array[array < 10e-10] = 0 # enough precision
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


def equalize(image):
    if isinstance(image, pil.Image):
        if image.mode in ("F"):
            return equalizearray(np.asarray(image))
        elif image.mode in ("RBGA"):
            image = image.convert("RBG")
        return ImageOps.equalize(image)
    else:
        return equalizearray(image)
