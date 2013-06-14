#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from bisect import insort

from blist import blist
import numpy as np
from numpy import NaN

from scipy.interpolate import griddata

from dft import get_sdct, get_idct, align_phase
from minimize import generic_minimizer
from image import get_subtract_paramns, normalize

pi = np.pi
tau = np.pi * 2 # two times sexier than pi


def get_aligned_phases(*phases):
    for phase in phases:
        uwphase = unwrap_wls(phase)
        if np.median(uwphase) > uwphase.mean():
            yield -phase
        else:
            yield phase


def compare_unwrappers(pea, unwrappers):
    unwrappeds = []
    for unwrapper in unwrappers:
        pea.unwrapper = unwrapper
        unwrappeds.append(normalize(pea.unwrapped_phase))
    return np.hstack(unwrappeds)


def wrapped_diff(phase, n=1, axis=-1, threshold=pi):
    diff = np.diff(phase, n, axis)
    diff[diff < -threshold] += 2 * threshold
    diff[diff > threshold] -= 2 * threshold
    return diff


def wrapped_gradient(phase):
    rows, cols = phase.shape
    dx, dy = np.gradient(phase)
    for diff in (dx, dy):
        diff[diff < -pi / 2] += pi
        diff[diff > pi / 2] -= pi

    gradient = dx + dy
    gradient[gradient < -pi / 4] = NaN
    gradient[gradient > pi / 4] = NaN

    gradient = fill_holes(gradient)

    return gradient


def fill_holes(array):
    array = array.copy()
    border = 3
    extended = np.hstack((array[:, -border:], array, array[:, :border]))
    extended = np.vstack((extended[-border:, :], extended, extended[:border, :]))

    nans = zip(*np.where(np.isnan(array)))
    print("Filling %d NaNs values with %d good data points." % 
        (len(nans), (array.size - len(nans))))
    for row, col in nans:
        context = extended[row: row + 2 * border + 1,
            col: col + 2 * border + 1]
        good_points = np.where(np.logical_not(np.isnan(context)))
        good_values = context[good_points]

        new_value = griddata(good_points, good_values, (border, border),
            method='cubic')
        array[row, col] = new_value
        extended[row + border, col + border] = new_value

    return array


def unwrap_wls(phase):
    """
    The fastest one but is innacurate.
    TODO: use lasso method
    """
    rows, cols = phase.shape

    wrowdiff = wrapped_diff(phase, 1, 0, pi)
    wrowdiff = np.concatenate((wrowdiff, np.zeros((1, cols))), 0)

    wcoldiff = wrapped_diff(phase, 1, 1, pi)
    wcoldiff = np.concatenate((wcoldiff, np.zeros((rows, 1))), 1)

    rhox = np.diff(np.concatenate((np.zeros((1, cols)),
        wrowdiff), 0), axis=0)
    rhoy =np.diff(np.concatenate((np.zeros((rows, 1)),
        wcoldiff), 1), axis=1)

    rho = rhox + rhoy
    dct = get_sdct(rho)

    col = np.mgrid[pi / cols:pi + pi / cols: pi / cols]
    row = np.mgrid[pi / rows:pi + pi / rows: pi / rows]
    cosines = 2 * (np.cos(row)[:, np.newaxis] + np.cos(col) - 2)

    try:
        phiinv = dct / cosines
    except:
        phiinv = dct / cosines[:-1, :-1]
    unwrapped = get_idct(phiinv)

    return unwrapped


def unwrap_wls2(phase):
    """
    The fastest one but is innacurate.
    TODO: use lasso method
    """
    rows, cols = phase.shape

    wrowdiff = wrapped_diff(phase, 1, 0, pi)
    wrowdiff = np.concatenate((wrowdiff, np.zeros((1, cols))), 0)

    wcoldiff = wrapped_diff(phase, 1, 1, pi)
    wcoldiff = np.concatenate((wcoldiff, np.zeros((rows, 1))), 1)

    rhox = np.diff(np.concatenate((np.zeros((1, cols)),
        wrowdiff), 0), axis=0)
    rhoy =np.diff(np.concatenate((np.zeros((rows, 1)),
        wcoldiff), 1), axis=1)

    rho = rhox + rhoy
    rho[rho < -pi / 4] = NaN
    rho[rho > pi / 4] = NaN

    rho = fill_holes(rho)
    dct = get_sdct(rho)

    col = np.mgrid[pi / cols:pi + pi / cols: pi / cols]
    row = np.mgrid[pi / rows:pi + pi / rows: pi / rows]
    cosines = 2 * (np.cos(row)[:, np.newaxis] + np.cos(col) - 2)

    try:
        phiinv = dct / cosines
    except:
        phiinv = dct / cosines[:-1, :-1]
    unwrapped = get_idct(phiinv)

    return unwrapped



def unwrap_multiphase(*phases):
    rows, cols = shape = phases[0].shape
    assert all((phase.shape == shape for phase in phases))

    rhos = []
    for phase in phases:
        wrowdiff = wrapped_diff(phase, 1, 0, pi)
        wrowdiff = np.concatenate((wrowdiff, np.zeros((1, cols))), 0)

        wcoldiff = wrapped_diff(phase, 1, 1, pi)
        wcoldiff = np.concatenate((wcoldiff, np.zeros((rows, 1))), 1)

        rhox = np.diff(np.concatenate((np.zeros((1, cols)),
            wrowdiff), 0), axis=0)
        rhoy =np.diff(np.concatenate((np.zeros((rows, 1)),
            wcoldiff), 1), axis=1)

        rhos.append(rhox + rhoy)

    rho = np.mean(rhos, 0)
    dct = get_sdct(rho)

    col = np.mgrid[pi / cols:pi + pi / cols: pi / cols]
    row = np.mgrid[pi / rows:pi + pi / rows: pi / rows]
    cosines = 2 * (np.cos(row)[:, np.newaxis] + np.cos(col) - 2)

    try:
        phiinv = dct / cosines
    except:
        phiinv = dct / cosines[:-1, :-1]
    unwrapped = get_idct(phiinv)

    return rho, unwrapped


def unwrap_qg(phase, quality_map):
    """
    Quality Guided Path Following unwrapping algoritm
    This algoritm uses the correlation array as quality map to guide the
    unwrapping path avoiding the tricky zones.

    Note: Correlation as also know as module image.

    Returns the unwrapped phase.
    """

    assert phase.shape == quality_map.shape
    phase = phase.copy()
    shape = phase.shape
    rows, cols = shape

    phase /= tau

    def get_neighbors(pos):
        row = pos / cols
        col = pos % cols
        if row > 0:
            yield pos - cols
        if row < (rows - 1):
            yield pos + cols
        if col > 0:
            yield pos - 1
        if col < (cols - 1):
            yield pos + 1

    phase = phase.ravel()
    adder = {}
    quality_map = quality_map.ravel()
    first_pixel = quality_map.argmax()
    border = blist()

    for pos in get_neighbors(first_pixel):
        adder[pos] = phase[first_pixel]
        insort(border, (quality_map[pos], pos))

    while border:
        quality, pixel = border.pop()
        phase[pixel] -= round(phase[pixel] - adder[pixel])

        for pos in get_neighbors(pixel):
            if pos not in adder:
                adder[pos] = phase[pixel]
                try:
                    insort(border, (quality_map[pos], pos))
                except IndexError:
                    print(quality_map.shape, pos)
                    raise

    phase = phase.reshape(shape) * tau
    return phase


def diff_match(phase1, phase2, threshold=np.pi/2.):
    """
    returns k that minimizes:

        var(diff(phase1) - diff(phase2 * k))

    """
    #TODO: re-implement minimization

    diffphase1 = wrapped_diff(phase1)
    diffphase2 = wrapped_diff(phase2)

    diffratio = diffphase1 / diffphase2
    diffratio[diffratio > threshold] = 1
    diffratio[diffratio < threshold ** -1] = 1

    best_k = diffratio.mean()

    print("var(bidiff(left) - bidiff(rigth * %f))" % best_k)
    return best_k


def phase_match(phase1, phase2):
    def diference(k):
        distance = ((phase1 - phase2 + k) ** 2).sum()
        return distance

    best_k = float(generic_minimizer(diference, 1))
    return best_k


def unwrap_phasediff(phase1, phase2):
    phase = phase1 - phase2
    phase[phase1 < phase2] += tau
    return phase


def unwrap_phasediff2(phase1, phase2):
    scale = diff_match(phase1, phase2)
    phase2 = phase2 * scale

    diff = np.zeros_like(phase1)
    diff[phase2 > phase1] += tau

    return phase2 + diff


def main():
    from pea import PEA
    from autopipe import showimage
    pea = PEA()
    pea.filename_holo = "../../../../documentos/carlos/enfused-sub/0427-0433-04X-568-c.tiff"
    wg = wrapped_gradient(pea.phase)
    showimage(wg)

    return 0


if __name__ == "__main__":
    exit(main())
