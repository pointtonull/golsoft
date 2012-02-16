#!/usr/bin/env python
#-*- coding: UTF-8 -*-

import os
import sys
import numpy as np


TKPIPE = os.name in ("nt")

if TKPIPE:
    import tkpipe
    import Image as pil
    TKPIPE = tkpipe.Tkpipe()
    sys.stdout = TKPIPE.default("green")
    sys.stderr = TKPIPE.default("red")
else:
    from scipy.misc import imshow


def showimage(image):
    if TKPIPE:
        if isinstance(image, np.ndarray):
            image = pil.fromarray(image)
        TKPIPE.writeimage(image)
        print("")
    else:
        imshow(image)

def color(message, color="blue"):
    if TKPIPE:
        TKPIPE.write(message, color)
    else:
        sys.stderr.write(message)

blue = lambda message:color(message, "blue")
red = lambda message:color(message, "red")
green = lambda message:color(message, "green")
