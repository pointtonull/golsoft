#!/usr/bin/env python
#-*- coding: UTF-8 -*-

import sys

from image import subtract, imread, imwrite
from os import path

def main():
    filenames = sys.argv[1:]
    if len(filenames) < 2:
        print("Se requieren al menos dos ficheros para operar la sustracción.")
        return 1

    files = [(filename, imread(filename)) for filename in filenames]
    name, left = files.pop(0)

    for right_name, right in files:
        left = subtract(left, right)
        filename_tail = path.split(right_name)[1]
        name += "-%s" % filename_tail

    print("Writting to %s" % name)
    imwrite(left, name)


if __name__ == "__main__":
    exit(main())
