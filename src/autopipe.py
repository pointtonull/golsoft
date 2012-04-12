#!/usr/bin/env python
#-*- coding: UTF-8 -*-


from StringIO import StringIO
import os
import sys


NP = False
try:
    import numpy as np
    NP = True
except:
    pass

ndarray = not NP or np.ndarray


if os.name in ("nt") or "TKPIPE" in os.environ:
    import tkpipe
    import Image as pil
    TKPIPE = tkpipe.Tkpipe()
    sys.stdout = TKPIPE.default("green")
    sys.stderr = TKPIPE.default("red")
else:
    from scipy.misc import imshow
    TKPIPE = False


def fig2raster(figure):
    """
    Convert a matplotlib to a raster PIL image
    """
    if hasattr(figure, "savefig"):
        fileo = StringIO()
        figure.savefig(fileo)
        fileo.seek(0)
        figure = pil.open(fileo)
    return figure


def showimage(*images):
    images = (fig2raster(image) for image in images)
    if TKPIPE:
        for image in images:
            if isinstance(image, ndarray):
                try:
                    image = pil.fromarray(image)
                except TypeError:
                    image = pil.fromarray(np.float64(image))
                except IndexError:
                    print image
                    raise
            TKPIPE.writeimage(image)
        print("")
    else:
        for image in images:
            imshow(image)
    return images


def color(message, color="blue"):
    if TKPIPE:
        TKPIPE.write(message, color)
    else:
        sys.stderr.write(message)
    return message


blue = lambda message:color(message, "blue")
red = lambda message:color(message, "red")
green = lambda message:color(message, "green")
