#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from bisect import insort

from blist import blist
import numpy as np

from dft import get_sdct, get_idct
from minimize import generic_minimizer

pi = np.pi
tau = np.pi * 2 # two times sexier than pi


def unwrap_wls(phase, quality_map=None):
    """
    Weighted least squares algoritm.
    The fastest one but is innacurate.
    TODO: implement weights!!
    TODO: try to use lasso method
    """
    rows, cols = phase.shape

    rowdiff = np.concatenate((np.diff(phase, 1, 0), np.zeros((1, cols))), 0)
    coldiff = np.concatenate((np.diff(phase, 1, 1), np.zeros((rows, 1))), 1)

    wrowdiff = np.mod(rowdiff + pi, tau) - pi
    wcoldiff = np.mod(coldiff + pi, tau) - pi

    rhox = np.diff(np.concatenate((np.zeros((1, cols)), wrowdiff), 0), axis=0)
    rhoy = np.diff(np.concatenate((np.zeros((rows, 1)), wcoldiff), 1), axis=1)

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


def get_bidiff(phase):
    if phase.ndim == 2:
        rows, cols = phase.shape
        rowdiff = np.concatenate((np.diff(phase, 1, 0), np.zeros((1, cols))),
            0)
        coldiff = np.concatenate((np.diff(phase, 1, 1), np.zeros((rows, 1))),
            1)

        wrowdiff = (rowdiff + pi) % tau - pi
        wcoldiff = (coldiff + pi) % tau - pi

        rhox = np.diff(np.concatenate((np.zeros((1, cols)), wrowdiff), 0),
            axis=0)
        rhoy =np.diff(np.concatenate((np.zeros((rows, 1)), wcoldiff), 1),
            axis=1)

        rho = rhox + rhoy
        return rho
    else:
        return np.diff(phase)


def get_aligned_phases(*phases):
    for phase in phases:
        uwphase = unwrap_wls(phase)
        if np.median(uwphase) > uwphase.mean():
            yield -phase
        else:
            yield phase


def unwrap_multiphase(*phases):
    rows, cols = shape = phases[0].shape
    assert all((phase.shape == shape for phase in phases))

    rhos = []
    for phase in phases:
        rho = get_bidiff(phase)
        rhos.append(rho)

    rho = np.median(rhos, 0)
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


def unwrap_multiphase2(*phases):
    rows, cols = shape = phases[0].shape
    assert all((phase.shape == shape for phase in phases))

    phases = [phase * diff_match(phases[0], phase)
        for phase in phases]

    rhos = []
    for phase in phases:
        rho = get_bidiff(phase)
        rhos.append(rho)

    rho = np.median(rhos, 0)
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
        row = pos % cols + 1
        col = pos / cols + 1
        if col > 1:
            yield pos - cols
        if col < cols:
            yield pos + cols
        if row > 1:
            yield pos - 1
        if row < rows:
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
                insort(border, (quality_map[pos], pos))

    phase = phase.reshape(shape) * tau
    return phase


def diff_match(phase1, phase2, threshold=np.pi/2.):
    """
    returns k that minimizes:

        var(diff(phase1) - diff(phase2 * k))
    """

    diffphase1 = get_bidiff(phase1)
    diffphase2 = get_bidiff(phase2)
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


def phasediff(phase1, phase2):
    phase = phase1 - phase2
    phase[phase1 < phase2] += tau
    return phase


def phasediff2(phase1, phase2):

    scale = diff_match(phase1, phase2)
    phase2 = phase2 * scale

    diff = np.zeros_like(phase1)
    diff[phase1 > phase2] += tau

    return phase2 + diff
