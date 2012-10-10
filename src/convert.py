#!/usr/bin/env python
#-*- coding: UTF-8 -*-

import sys

from image import imread, imwrite


def main():
    filenames = [filename for filename in sys.argv[1:]]
    assert len(filenames) == 2
    fromfile, tofile = filenames
    print("Converting '%s' to '%s'" % (fromfile, tofile))
    image = imread(fromfile)
    imwrite(image, tofile)


if __name__ == "__main__":
    exit(main())
