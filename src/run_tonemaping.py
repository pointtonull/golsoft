"""
A simple interactive amoeba demo
"""

from ibisect import Amoeba
from image import imread, imwrite, equalize, normalize
from scipy.ndimage.filters import median_filter, gaussian_filter
from scipy import misc
import numpy as np
import sys

from autopipe import showimage

def main():
    """
    Le main rutine
    """
    #TODO: se puede automatizar los mejores paramétros para quitar las pertubaciones del fondo
    images = [(filename, imread(filename, False))
        for filename in sys.argv[1:]]
    if not images:
        lena = misc.lena()
        images = [("lena", lena)]

    def sigmoid(x):
        result = np.abs(1 / (1 + np.exp(-x)))
        return result

    for filename, image in images:
        print("Original %s:" % filename)

        def processor(t_sigma, t_level, equalize_level):
            t_sigma = sigmoid(t_sigma) * 20 + 1
            t_level = sigmoid(t_level)
            equalize_level = sigmoid(equalize_level)
            channels = []
            for ndim in range(3):
                channel = image[:, :, ndim]
                local_context = gaussian_filter(channel, t_sigma)
#                local_context = median_filter(channel, t_sigma)
                tonemapped = channel - local_context * t_level
                tonemapped = tonemapped.astype(int)
                tonemapped /= 1 - t_level
                equalized = equalize(tonemapped) * equalize_level
                equalized += tonemapped * (1 - equalize_level)
                channels.append(equalized)
            final = np.array(channels).swapaxes(0, 1).swapaxes(1, 2)
            final = normalize(final).astype(np.uint8)
            return final

        amoeba = Amoeba(processor)
        amoeba.iterate(distance=1)
        print(amoeba.points)

if __name__ == "__main__":
    exit(main())
