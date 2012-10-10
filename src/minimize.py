#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from scipy import optimize

def generic_minimizer(fitness_func, initial_guess, optimizers=None):
    """
    A common interface to various minimization algorithms
    """

    if optimizers == None:
        optimizers = [
            optimize.fmin, # 66
            optimize.fmin_powell,
        ]

    best_result = None
    for optimizer in optimizers:
        xend = optimizer(fitness_func, initial_guess, disp=False)
        last_result = fitness_func(xend)
        if best_result is None or last_result < best_result:
            best_guess = xend
            best_result = last_result

    return best_guess

