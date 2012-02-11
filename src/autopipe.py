#!/usr/bin/env python
#-*- coding: UTF-8 -*-

import os
import sys


TKPIPE = os.name in ("nt")

if TKPIPE:
    import tkpipe
    TKPIPE = tkpipe.Tkpipe()
    sys.stdout = pipe.default("green")
    sys.stderr = pipe.default("red")


def showimage(image):
    if TKPIPE:
        TKPIPE.writeimage(image)
    else:
        image.show()

def color(message, color="blue"):
    if TKPIPE:
        TKPIPE.write(message, color)
    else:
        sys.stderr.write(message)

blue = lambda message:color(message, "blue")
red = lambda message:color(message, "red")
green = lambda message:color(message, "green")
