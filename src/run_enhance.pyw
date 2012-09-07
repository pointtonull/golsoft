#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from autopipe import showimage
from image import imread, derotate, get_centered, phase_denoise, normalize
from scipy.misc import imresize, imrotate, lena
from itertools import product
from random import sample
import sys
import numpy as np

tau = np.pi * 2


def main():
    images = [(filename, get_centered(imread(filename)))
        for filename in sys.argv[1:]]

    if not images:
        image = imresize(lena(), .5)
        images = [("lena", image)]

    for name, image in images:
        print(name)
        showimage(image)
        phase = (image / 10.) % tau
        showimage(normalize(phase))
        print phase.ptp()
        denoised = phase_denoise(phase)
        showimage(normalize(denoised))
        print denoised.ptp()

    return 0

if __name__ == "__main__":
    exit(main())
