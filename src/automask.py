#!/usr/bin/env python
#-*- coding: UTF-8 -*-

"""
Simple module to apply automatic masks to select the correct orders on a
spectrum image.
View run_automask to see a complete example.
"""

from scipy.ndimage import geometric_transform
from scipy.ndimage import maximum_filter1d, rotate
import numpy as np

tau = np.pi * 2
HOLE_LEN = 4


def get_localmaxs(array, count=3, window=25):
    """
    Individualize the locals maxs relatives to a window of the given length.
    Returns theirs coodenates.
    """
    maxs = maximum_filter1d(array, window)
    maxs[maxs > array] = array.min()
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


def get_circles(array, count=3, window=30):
    """
    Identify the center and radius for each non 0 order.
    Returns: value, center, radius
    Value: the value at the center cell
    Center: row, col position
    Radius: ¬¬
    """
    #FIXME: too large and bored
    rowsums = array.sum(1)
    rowpeaks = get_localmaxs(rowsums, count, window)
    rowsstd = np.std([rowsums[pos] for pos in rowpeaks])

    colsums = array.sum(0)
    colpeaks = get_localmaxs(colsums, count, window)
    colsstd = np.std([colsums[pos] for pos in colpeaks])

    if colsstd < rowsstd:

        centers = []
        for colpeak in colpeaks:
            rowsums = array[:, colpeak - 5:colpeak + 5].sum(1)
            rowpeak = get_localmaxs(rowsums, 1)[0]
            centers.append((int(rowpeak), int(colpeak)))

    else:

        centers = []
        for rowpeak in rowpeaks:
            colsums = array[rowpeak - 5:rowpeak + 5].sum(0)
            colpeak = get_localmaxs(colsums, 1)[0]
            centers.append((int(rowpeak), int(colpeak)))

    circles = [] #value, center, radio
    valley = get_wide_segment(array, centers[0], centers[1])
    peaks = get_localmaxs(-valley, 1)[0]
    radius = abs(int(round(peaks)))
    circles.append((array[centers[0]], centers[0], radius))

    valley = get_wide_segment(array, centers[1], centers[2])
    peaks = get_localmaxs(-valley, 1)[0]
    radius = abs(int(round(centers[2][0] - centers[1][0] - peaks)))
    circles.append((array[centers[2]], centers[2], radius))

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


def get_centered(array, center):
    """
    Shift the given array to make the given point be the new center.
    """
    #TODO: can improve speed using ndimage.shift
    rows, cols = array.shape
    rowcc = rows / 2.
    colcc = cols / 2
    rowc, colc = center
    drows = rowc - rowcc 
    dcols = colc - colcc 
    
    def out2in((outrow, outcol)):
        return outrow + drows, outcol + dcols

    centered = geometric_transform(array, out2in, array.shape, 
        mode="wrap")
    return centered


def get_mask(shape, window, center=None):
    """
    Draw a mask of the given shape and window.
    If center is not given it will be the shape center.
    """
    array = np.zeros(shape)
    window2d = radial_extrusion(window)
    if center is None:
        center = shape[0] / 2., shape[1] / 2.
    crow, ccol = center
    radius = window.shape[0] / 2.
    top, bottom = (crow - radius), (crow + radius)
    left, right = (ccol - radius), (ccol + radius)
    array[top:bottom, left:right] = window2d
    return array


def get_holed_window(winfunc, length, holelen=4):
    """
    Create a window with a centered hole.
    """
    window = winfunc(length)
    center = length / 2.
    holelen -= holelen % 2 #FIXME: only evens for now
    hole = np.ones(holelen) - winfunc(holelen)
    window[center - holelen / 2:center + holelen / 2] *= hole
    return window



def main():
    print("Import then exist.")


if __name__ == "__main__":
    exit(main())
