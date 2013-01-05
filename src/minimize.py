#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from scipy import optimize, stats
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
    return abs(((array1 - array2) ** 2)).sum()


def get_paraboloid(x, y, a0, b0, a1, b1, c=0):
    """ a0 * (x - b0) ** 2 + a1 * (y - b1) ** 2 + c """
    return a0 * (x - b0) ** 2 + a1 * (y - b1) ** 2 + c


def get_plane(linspace, a, b):
    """ linspace * a + b """
    return linspace * a + b


def get_fitted_paraboloid(data):
    xs, ys = data.shape
    x, y = np.mgrid[:xs, :ys]

    def fitness((a0, b0, a1, b1, c)):
        error = squared_error(data, get_paraboloid(x, y, a0, b0, a1, b1, c))
        return error

    params = generic_minimizer(fitness, [1] * 5)
    return get_paraboloid(x, y, *params)


def wrapped_gradient(phase):
    rows, cols = phase.shape
    dx, dy = np.gradient(phase)
    for diff in (dx, dy):
        diff[diff < -pi / 2] += pi
        diff[diff > pi / 2] -= pi

    gradient = dx + 1j * dy

    return gradient


def get_fitted_paraboloid2(data):
    """
    Similar to get_fitted_paraboloid but uses the complex gradient to determine
    the fittness. This method allow us to correct a wrapped phase.
    """
    xs, ys = data.shape
    x = np.mgrid[:xs]
    y = np.mgrid[:ys]

    diff_outline_x = np.gradient(data.mean(1))
    diff_outline_y = np.gradient(data.mean(0))

    dax, dbx, r_value, p_value, std_err = stats.linregress(x, diff_outline_x)
    day, dby, r_value, p_value, std_err = stats.linregress(y, diff_outline_y)

    ax = dax / 2 
    bx = - dbx / dax
    ay = day / 2 
    by = - dby / day
    x, y = np.mgrid[:xs, :ys]
    return get_paraboloid(x, y, ax, bx, ay, by)


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
