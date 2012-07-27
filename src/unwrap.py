#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from bisect import insort

from blist import blist
import image
import numpy as np

from dft import get_sdct, get_idct

pi = np.pi
tau = 2 * pi


def unwrap_wls(wrapped):
    """
    Weighted least squares algoritm.
    The fastest one but extremly innacurate.
    TODO: implement weights!!
    TODO: try to use lasso method
    """
    rows, cols = wrapped.shape

    rowdiff = np.concatenate((np.diff(wrapped, 1, 0), np.zeros((1, cols))), 0)
    coldiff = np.concatenate((np.diff(wrapped, 1, 1), np.zeros((rows, 1))), 1)

    wrowdiff = np.mod(rowdiff + pi, tau) - pi
    wcoldiff = np.mod(coldiff + pi, tau) - pi

    rhox = np.diff(np.concatenate((np.zeros((1, rows)), wrowdiff), 0), axis=0)
    rhoy = np.diff(np.concatenate((np.zeros((cols, 1)), wcoldiff), 1), axis=1)

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


def unwrap_qg(phase, quality_map):
    """
    Quality Guided Path Following unwrapping algoritm
    This algoritm uses the correlation array as quality map to guide the
    unwrapping path avoiding the tricky zones.

    Note: Correlation as also know as module image.

    Returns the unwrapped phase.
    """
    
    assert phase.shape == quality_map.shape
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


"""
Implementarlo ahora se sale del cronograma pero creo que tengo una idea que
captura la misma información que multi-path quality guided unwrap y es muy
vectorizable.

    1. Crear un mapa de pendientes a partir del mapa de calidad. Este mapa de
       pendientes codifica de dirección en que se debe derivar e integrar. En cada
       pixel codifica guarda un número del [1, 4] donde 1 es arriba, 2, derecha,
       ect. La dirección a escoger debe ser siempre la que cree la de mejor
       calidad (el vecino con mayor valor en el mapa de calidad).
    2. Crear la derivada discreta siguiendo el mapa de pendientes.
    3. Modificar la derivada para que todos los números de acerquen lo más posible
       a 0 eliminando las congruensias mayores. matriz =% pi ¿?
    4. Integrar 3. siguiendo el mapa de pendientes.

Lo truculento está en la integración. Todos los pasos previos son sencillos y
relativamente faciles de programar. Al integrar de puede partir del máximo
absoluto, asignarle 0. Luego integrar sus vecinos en orden de peor a mejor
calidad (para evitar que los pixeles con mala calidad sobreescriban a los que
tienen buena calidad). Repetir hasta haber integrado todo el mapa. Las áreas del
mapa con mala calidad se aislan solas rompiendo el avance uniforme pero caen al ser
rodeadas.
"""
