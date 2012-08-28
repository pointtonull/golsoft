#!/usr/bin/env python
#-*- coding: UTF-8 -*-

"""
Simple module to apply automatic masks to select the correct orders on a
spectrum image.
View run_automask to see a complete example.
"""

#from autopipe import showimage
#from image import normalize
from collections import defaultdict
from image import get_centered, get_intensity
from scipy.ndimage import geometric_transform
from scipy.ndimage import maximum_filter1d, rotate
import matplotlib.pyplot as plt
import numpy as np

tau = np.pi * 2
VERBOSE = 0


def graph(*arrays):
    if VERBOSE:
        for array in arrays:
            plt.plot(array)
        plt.show()
        plt.close()


def get_localmaxs(array, count=3, window=25):
    """
    array must be one-dimensional
    Individualize the locals maxs relatives to a window of the given length.
    Returns theirs coodenates.
    """
    maxs = maximum_filter1d(array, window)
    maxs[maxs > array] = array.min()
    graph(array, maxs)
    top_positions = maxs.argsort()[::-1][:count]
    return sorted(top_positions)


def get_wide_segment(array, startpoint, endpoint):
    """
    Returns an array with the values inside of the 2x1 rectangle so that
    the centers of the shorter sides are at the start and end point. 
    """
    row0, col0 = startpoint
    row1, col1 = endpoint
    drow = row1 - row0
    dcol = col1 - col0
    hyp = (drow ** 2 + dcol ** 2) ** .5

    center = (row0 + row1) / 2., (col0 + col1) / 2.
    centered = get_centered(array, center)
    rotated = rotate(centered, np.arctan2(drow, dcol) * 360 / tau, mode='wrap')
    rrows, rcols = rotated.shape
    rrowc = rrows / 2.
    rcolc = rcols / 2.
    col0 = rcolc - hyp / 2.
    col1 = rcolc + hyp / 2.
    croped = rotated[rrowc - hyp / 2: rrowc + hyp / 2, col0:col1]
    sumed = croped.sum(0)
    return sumed


def get_circles(array, count=3, window=25):
    """
    Identify the center and radius for each local max
    Returns value, center, radius
        Value: the value at the center cell
        Center: row, col position
        Radius: radius
    """

    rowsums = array.sum(1)
    rowpeaks = get_localmaxs(rowsums, count, window)
    rowsstd = np.std([rowsums[pos] for pos in rowpeaks])

    colsums = array.sum(0)
    colpeaks = get_localmaxs(colsums, count, window)
    colsstd = np.std([colsums[pos] for pos in colpeaks])

    centers = []
    if colsstd < rowsstd:
        for colpeak in colpeaks:
            rowsums = array[:, colpeak - 5:colpeak + 5].sum(1)
            rowpeak = get_localmaxs(rowsums, 1)[0]
            centers.append((int(rowpeak), int(colpeak)))
    else:
        for rowpeak in rowpeaks:
            colsums = array[rowpeak - 5:rowpeak + 5].sum(0)
            colpeak = get_localmaxs(colsums, 1)[0]
            centers.append((int(rowpeak), int(colpeak)))

    radius = defaultdict(float)
    for center0, center1 in zip(centers, centers[1:]):
        valley = get_wide_segment(array, center0, center1)
        peaks = get_localmaxs(-valley, 1)[0]
        radius0 = abs(int(round(peaks)))
        radius1 = valley.shape[0] - radius0
        radius[center0] = max(radius[center0], radius0)
        radius[center1] = max(radius[center1], radius1)

    circles = [(array[center], center, radius[center]) for center in centers]
    return circles


def radial_extrusion(array, center=None):
    """
    Create a radial extrusion from one-dimensional array passed as an argument.
    If not specified the center of rotation will be set to the center of array.
    """
    inshape = array.shape
    assert len(inshape) == 1
    center = center or inshape[0] / 2.
    outxs = max(inshape[0] - center, center) * 2
    outys = outxs

    def out2in((outx, outy)):
        rho = ((outx-center) ** 2 + (outy-center) ** 2) ** .5
        return (center + rho, )

    extrusion = geometric_transform(array, out2in, (outxs, outys))
    return extrusion


def get_mask(shape, window, center=None):
    """
    Draw a mask of the given shape and window.
    If center is not given it will be the shape center.
    """

    array = np.zeros(shape)
    window2d = radial_extrusion(window)

    array_center = shape[0] / 2., shape[1] / 2.
    c_row, c_col = array_center
    radious = window.shape[0] / 2.
    top, bottom = (c_row - radious), (c_row + radious)
    left, right = (c_col - radious), (c_col + radious)

    try:
        array[top:bottom, left:right] = window2d
    except ValueError:
        print center
        print array[top:bottom, left:right].shape, window2d.shape
        print "array[%d:%d, %d:%d] = [](%s)" % (top, bottom, left, right, 
            window2d.shape)
        raise

    if center:
        array = get_centered(array, center, reverse=True)

    return array


def get_holed_window(winfunc, length, holelen=0):
    """
    Create a window with a centered hole.
    """
    if holelen > length:
        holelen = 0
    length = int(round(length))
    window = winfunc(length)
    center = length / 2.
    if holelen:
        holelen -= holelen % 2 #FIXME: only evens for now
        hole = 1 - winfunc(holelen)
        try:
            window[center - holelen / 2:center + holelen / 2] *= hole
        except ValueError:
            print(length, holelen)
            raise
    return window


def get_auto_mask(spectrum, softness=1, radious_scale=1, zero_scale=1,
        cuttop=0):
    """
    Try to filter spurious data out.
    """
    shape = spectrum.shape
    intensity = get_intensity(spectrum)

    circles = sorted(get_circles(intensity, 3, 50))
    virtual_order, real_order, zero_order = circles
    peak_height, peak_center, peak_radious = real_order

    peak_radious = min([(abs(shape[0] / 3.5 - peak[2]), peak[2])
        for peak in circles])[1]

    windowmaker = lambda x: np.kaiser(x, softness)
    window = get_holed_window(windowmaker, peak_radious * radious_scale)
    mask = get_mask(shape, window, peak_center)

    zerowindow = get_holed_window(windowmaker, peak_radious * zero_scale)
    zeromask = 1 - get_mask(shape, zerowindow, zero_order[1])
    mask *= zeromask

    masked_intensity = mask * intensity

    cutoff = masked_intensity > (masked_intensity.max()
        - masked_intensity.ptp() * cuttop)
    mask[cutoff] = 0
    masked = mask * spectrum

    centered = get_centered(intensity, peak_center)
    masked = get_centered(masked, peak_center)
    mask = get_centered(mask, peak_center)

    return mask, masked, centered



def main():
    import sys
    from scipy import misc
    softness = 2
    r_scale = 3
    windowmaker = lambda x: np.kaiser(x, softness)
    for infilename in sys.argv[1:]:
        print("File: %s" % infilename)
        array = misc.imread(infilename)
        circles = get_circles(array, 5)
        l_pnoise, l_order, c_order, r_order, r_pnoise = circles

        window = get_holed_window(windowmaker, l_order[2] * r_scale, 10)
        mask = get_mask(array.shape, window, l_order[1])
        window = windowmaker(l_pnoise[2] * r_scale)
        mask *= 1 - get_mask(array.shape, window, l_pnoise[1])
        window = windowmaker(c_order[2] * r_scale)
        mask *= 1 - get_mask(array.shape, window, c_order[1])
#        showimage(mask * 255)
#        showimage(logscale(mask * array))

        window = get_holed_window(windowmaker, r_order[2] * r_scale, 10)
        mask = get_mask(array.shape, window, r_order[1])
        window = windowmaker(r_pnoise[2] * r_scale)
        mask *= 1 - get_mask(array.shape, window, r_pnoise[1])
        window = windowmaker(c_order[2] * r_scale)
        mask *= 1 - get_mask(array.shape, window, c_order[1])

#        showimage(mask * 255)
#        showimage(logscale(mask * array))


if __name__ == "__main__":
    exit(main())
