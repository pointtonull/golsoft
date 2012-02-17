#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from autopipe import showimage, red, blue
from enhance import equalize, autocontrast, logscale
from itertools import combinations
from scipy import misc
from scipy.ndimage import geometric_transform, gaussian_filter
from scipy.ndimage import maximum_filter1d, minimum_filter1d, rotate
import Image as pil
import numpy as np
#import pylab
import scipy
import sys

tau = np.pi * 2
HOLE_LEN = 4


def get_localmaxs(array, count=3, window=25):
    # just a concept, must be improved
    maxs = maximum_filter1d(array, window)
    maxs[maxs > array] = array.min()
    top_positions = maxs.argsort()[::-1][:count]
    return sorted(top_positions)


def get_peaks1d(array, count=1, window=30):
    localmaxs = get_localmaxs(array, count, window)
    return localmaxs


def graf(*data):
    for datum in data:
        pylab.plot(datum)
    return pylab.show()


def get_segment(array, startpoint, endpoint):
    startingrow, startingcol = startpoint
    endingrow, endingcol = endpoint
    rowtravel = endingrow - startingrow
    coltravel = endingcol - startingcol
    totaltravel = (rowtravel ** 2 + coltravel ** 2) ** .5

    def out2in(hipotenusa):
        rowfrom = startingrow + rowtravel * hipotenusa / totaltravel
        colfrom = startingcol + coltravel * hipotenusa / totaltravel
        return rowfrom, colfrom

    segment = geometric_transform(array, out2in, (totaltravel, ))
    return segment


def get_wide_segment(array, startpoint, endpoint):
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
    rowsums = array.sum(1)
    rowpeaks = get_peaks1d(rowsums, count, window)
    rowsstd = np.std([rowsums[pos] for pos in rowpeaks])

    colsums = array.sum(0)
    colpeaks = get_peaks1d(colsums, count, window)
    colsstd = np.std([colsums[pos] for pos in colpeaks])

    if colsstd < rowsstd:

        centers = []
        for colpeak in colpeaks:
            rowsums = array[:, colpeak - 5:colpeak + 5].sum(1)
            rowpeak = get_peaks1d(rowsums, 1)[0]
            centers.append((int(rowpeak), int(colpeak)))

    else:

        centers = []
        for rowpeak in rowpeaks:
            colsums = array[rowpeak - 5:rowpeak + 5].sum(0)
            colpeak = get_peaks1d(colsums, 1)[0]
            centers.append((int(rowpeak), int(colpeak)))

    print centers
    circles = [] #value, center, radio
    valley = get_wide_segment(array, centers[0], centers[1])
    peaks = get_peaks1d(-valley, 1)[0]
    radius = abs(int(round(peaks)))
    circles.append((array[centers[0]], centers[0], radius))

    valley = get_wide_segment(array, centers[1], centers[2])
    peaks = get_peaks1d(-valley, 1)[0]
    radius = abs(int(round(centers[2][0] - centers[1][0] - peaks)))
    circles.append((array[centers[2]], centers[2], radius))

    return circles


def draw_circle(canvas_shape, center, radius, fill=1):
    xc, yc = center
    r = radius
    array = np.ndarray(canvas_shape)
    Yiter, Xiter = np.ogrid[0:canvas_shape[0] , 0:canvas_shape[1]]
    mask = (Xiter - xc) ** 2 + (Yiter - yc) ** 2 <= r ** 2
    array[mask] = fill
    return array


def radial_extrusion(array, center=None):
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


def kaiser14(length):
    return np.kaiser(length, 14)


def get_holed_mask(winfunc, length, holelen=4):
    window = winfunc(length)
    center = length / 2.
    holelen -= holelen % 2 #only odds for now
    hole = np.ones(holelen) - winfunc(holelen)
    window[center - holelen / 2:center + holelen / 2] *= hole
    return window



def main():

    windows = {
        "Rectangular": np.ones,
        "Hanning": np.hanning,
        "Bartlett": np.hamming,
        "Blackman": np.bartlett,
        "Hamming": np.blackman,
        "Kaiser14": kaiser14,
    }

    if len(sys.argv) < 2:

        for w_name in windows:

            print("%s window:" % w_name)
            window = windows[w_name](500)
            window2d = get_mask((500, 500), window)
            showimage(window2d * 255)

    else:

        for infilename in sys.argv[1:]:
            print("File: %s" % infilename)
            array = misc.imread(infilename)
            orows, ocols = array.shape
            if "." in infilename:
                infilename = "".join(infilename.split(".")[:-1])
            circles = sorted((get_circles(array)))
            for value, (crow, ccol), radius in circles:
                print("  Center %d, %d, radius %d:" % (ccol, crow, radius))
                top, bottom = crow - radius, crow + radius
                left, right = ccol - radius, ccol + radius
                helprect = "X%d-%d Y%d-%d" % (left, right, top, bottom)
                print("    %s" % helprect)
                for w_name, w_function in windows.iteritems():
                    plusholes = {
                        "o" + w_name: lambda length:get_holed_mask(w_function,
                            length, HOLE_LEN),
                        w_name: w_function,
                    }
                    for w_name, w_func in plusholes.iteritems():
                        filename = "-".join((infilename, w_name, helprect))
                        filename += ".tiff"
                        print("      %s" % w_name)
                        window = w_func(radius * 2)
                        mask = get_mask((512, 512), window, (crow, ccol))
                        image = pil.fromarray(np.float32(mask))
                        image.save(filename)
                        showimage(mask * 255)



if __name__ == "__main__":
    exit(main())
