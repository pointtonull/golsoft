#!/usr/bin/env python

import cv2 as cv
#import mdp
import numpy as np
import unwrap

from autopipe import showimage
from image import normalize
from numpy import hstack, vstack
from pea import PEA


def main():
    pea = PEA()
    pea.filename_holo = "1206-h-det.png"
    pea.filename_obj  = "1206-o-det.png"
    pea.filename_ref  = "1206-r-det.png"
    pea.unwrapper = unwrap.unwrap_qg
    phase = normalize(pea.unwrapped_phase)
    module = pea.image_obj

    matrix = vstack((phase.flatten(), module.flatten()))

    print("Computing PCA")
    mean, eigenvectors = cv.PCACompute(matrix, np.mean(matrix, axis=0).reshape(1,-1))
    print(type(eigenvectors))
    print(dir(eigenvectors))
    print("Projecting data")
    pca = cv.PCAProject(matrix, mean, eigenvectors)
    print(dir(pca))
    return
    backimage = cv.PCABackProject(pca, mean, eigenvectors)
    showimage(backimage.reshape(phase.shape))


    return 0

if __name__ == "__main__":
    exit(main())
