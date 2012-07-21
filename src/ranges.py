#!/usr/bin/env python
#-*- coding: UTF-8 -*-


def frange(mean, radius, amount):
    """
    Helper for binary search of floats values
    """
    if amount == 1:
        start = mean
        step = 0.
    else:
        start = mean - radius
        step = (radius * 2) / (amount - 1.)
    for stepn in xrange(amount):
        yield start + stepn * step



def main():
    for i, n in enumerate(frange(0, 5, 20)):
        print i, n

if __name__ == "__main__":
    exit(main())
