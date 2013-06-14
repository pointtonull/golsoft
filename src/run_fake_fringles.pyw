#!/usr/bin/env python
#-*- coding: UTF-8 -*-

import sys

from numpy import hstack, pi, fabs, sum
from pea import get_initial_phase, get_fringes
from scipy import misc, ndimage, optimize
import numpy as np

from automask import get_mask
from minimize import generic_minimizer
from autopipe import showimage
from image import imread, normalize

tau = 2 * pi


def get_similarity_map(left, rigth, blur=10):
    overlap = - fabs(left - rigth)
    overlap = ndimage.filters.gaussian_filter(overlap, blur)
    overlap -= overlap.min()
    overlap /= overlap.max()
    return overlap


def get_patch(shape, row, col, height, diameter):
    if row < 0:
        row = 0
    elif row >= shape[0]:
        row = shape[0] - 1

    if col < 0:
        col = 0
    elif col >= shape[1]:
        col = shape[1] - 1

    min_side = min(shape)
    if diameter > min_side:
        diameter = min_side
    window = np.hanning(round(diameter)) * height
    patch = get_mask(shape, window, (row, col))
    return patch


def main():
    images = [(filename, imread(filename, True)) for filename in sys.argv[1:]]
    if len(images) < 2:
        if not images:
            lena = misc.imresize(misc.lena(), .5)
            images = [("lena", lena)]

    for filename, hologram in images:
        print("Original hologram: %s" % filename)
        hologram = normalize(hologram)
        print("shape: %d, %d" % hologram.shape)
        # procesar para regularizar iluminación y quitar algo de ruido
        hologram -= hologram.min()
        hologram /= hologram.max()
        hologram = hologram * 2 - 1 # a la imagen de seno

        phase = get_initial_phase(hologram)

        SIMILARITY_THRESHOLD = 0.9

        for step in xrange(500):
            fringes = get_fringes(phase)

            similarity_map = get_similarity_map(hologram, fringes)

            if not step % 5:
                showimage(hstack((hologram, fringes, (2 * similarity_map) - 1)))

            similarity_map[similarity_map > SIMILARITY_THRESHOLD] = 0

            if not np.any(similarity_map):
                break # ya no quedan zonas a mejorar
            row, col = np.unravel_index(similarity_map.argmax(), phase.shape)

            def fitness_func((height, diameter)):
                print("row: %.2f, col: %.2f, height: %.2f, diameter: %.2f"
                    % (row, col, height, diameter))
                patch = get_patch(phase.shape, row, col, height, diameter)
                new_phase = phase + patch
                new_fringes = get_fringes(new_phase)
                score = sum(np.fabs(hologram - new_fringes))
                print("score: %.2f" % score)
                return score
            initial_guess = np.array([0, 30], float)
            print(initial_guess)
            best_change = generic_minimizer(fitness_func, initial_guess, [optimize.fmin])
            print(best_change)
            patch = get_patch(phase.shape, row, col, *best_change)
            phase += patch

        showimage(phase)

    return 0


if __name__ == "__main__":
    exit(main())
