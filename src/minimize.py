#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from scipy import optimize
from numpy import pi

import numpy as np

tau = 2 * pi

def generic_minimizer(fitness_func, initial_guess, optimizers=None):
    """
    A common interface to various minimization algorithms
    """
    
    if optimizers == None:
        optimizers = [
            optimize.fmin, # 66
            optimize.fmin_powell,
#            optimize.leastsq,
        ]

    best_result = None
    for optimizer in optimizers:
        xend = optimizer(fitness_func, initial_guess, disp=False)
        last_result = fitness_func(xend)
        if best_result is None or last_result < best_result:
            best_guess = xend
            best_result = last_result

    return best_guess


def squared_error(array1, array2):
    return ((array1 - array2) ** 2).sum()


def get_paraboloid(x, y, a0, b0, a1, b1, c):
    return (a0 * x + b0) ** 2 + (a1 * y + b1) ** 2 + c


def get_fitted_paraboloid(data):
    xs, ys = data.shape
    x, y = np.mgrid[:xs, :ys]

    def fitness((a0, b0, a1, b1, c)):
        error = squared_error(data, get_paraboloid(x, y, a0, b0, a1, b1, c))
        return error

    params = generic_minimizer(fitness, [1] * 5)
    return get_paraboloid(x, y, *params)


def main():
    from scipy.misc import lena
    from autopipe import showimage
    x, y = np.mgrid[:100, :100]
    eye = lena()[200:300, 200:300]
    data = get_paraboloid(x, y, 1, 2, 3, 4, 5)
    noisy = data + eye
    fitted = get_fitted_paraboloid(noisy)
    showimage(eye, noisy, fitted, noisy - fitted)
    return 0


if __name__ == "__main__":
    exit(main())
