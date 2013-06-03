#!/usr/bin/env python
#-*- coding: UTF-8 -*-

import sys

from numpy import mean

from autopipe import showimage
from image import normalize, imread
from pea import PEA
from unwrap import unwrap_qg


def main():
    """
    Le main rutine
    """
    filenames = sys.argv[1:]
    pea = PEA()
    if len(filenames) == 9:
        images = [imread(filename) for filename in filenames]

        holo = mean(images[:3], 0)
        obj = mean(images[3:6], 0)
        ref = mean(images[6:], 0)

        pea.image_holo = holo
        pea.image_obj = obj
        pea.image_ref = ref
    elif len(filenames) == 3:
        pea.filename_holo = filenames[0]
        pea.filename_obj = filenames[1]
        pea.filename_ref = filenames[2]
    elif len(filenames) == 1:
        pea.filename_holo = filenames[0]
    else:
        print("la cantidad de argumentos debe ser 1, 3 o 9")
        return 1

    showimage(pea.image)

    pea.propagate = False
    showimage(normalize(pea.phase))
    showimage(normalize(pea.module))
    showimage(normalize(pea.unwrapped_phase))
    pea.unwrapper = unwrap_qg
    showimage(normalize(pea.unwrapped_phase))


    return 0


if __name__ == "__main__":
    exit(main())
