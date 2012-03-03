#!/usr/bin/env python
#-*- coding: UTF-8 -*-

import time
import sys

INICIO = time.time()
VERBOSE = 1

class Verbose:
    def __init__(self, verbosity, prefix="", ident=True):
        self.verbosity = False if verbosity < 0 else True
        self.prefix = prefix
        self.ident = ident


    def __call__(self, *args):
        if self.verbosity:
            message = " ".join((unicode(e) for e in args))
            sys.stderr.write("%s%s%s\n" % (" " * self.get_depth(), self.prefix,
                message))


    def get_depth(self):
        if not self.ident:
            return 0
        else:
            def exist_frame(n):
                try:
                    if sys._getframe(n):
                        return True
                except ValueError:
                    return False

            now = 0
            maxn = 1
            minn = 0

            while exist_frame(maxn):
                minn = maxn
                maxn *= 2

            middle = (minn + maxn) / 2
          
            while minn < middle:
                if exist_frame(middle):
                    minn = middle
                else:
                    maxn = middle

                middle = (minn + maxn) / 2
          
            return max(minn - 3, 0) #4 == len(main, Verbose, get_depth)


error = Verbose(VERBOSE + 2, "E: ")
warning = Verbose(VERBOSE + 1, "W: ")
info = Verbose(VERBOSE + 0)
moreinfo = Verbose(VERBOSE -1)
debug = Verbose(VERBOSE - 2, "D: ")
