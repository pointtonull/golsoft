#!/usr/bin/env python
#-*- coding: UTF-8 -*-

import sys

from autopipe import showimage
from image import normalize
from pea import PEA


def main():
    """
    Le main rutine
    """
    filenames = sys.argv[1:]

    for filename in filenames:
        print(filename)
        pea = PEA(filename)
        showimage(normalize(pea.image))

        print("Using manual focus")
        pea.propagate = False
        showimage(normalize(pea.phase))
        showimage(normalize(pea.module))
        showimage(normalize(pea.unwrapped_phase))

    return 0


if __name__ == "__main__":
    exit(main())
