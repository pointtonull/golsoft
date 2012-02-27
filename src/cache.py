#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from functools import wraps
import inspect
import os
import pickle
import sys
import time
import debug
import hashlib

debug.VERBOSE = 0
debug = debug.debug


def _hash_arg(args):
    try:
        hash(args)
        return args
    except TypeError:
        try:
            ahash = hashlib.sha1(args).hexdigest()
            return ahash
        except TypeError:
            values = tuple((_hash_arg(arg) for arg in args))
            return values


class Cache:
    def __init__(self, filename=None, deadline=100*86400, flush_frequency=0):
        self.count = 0
        self.deadline = deadline
        self.filename = filename
        self.flush_frequency = flush_frequency

        if filename:
            try:
                self.cache = pickle.load(open(self.filename))
            except IOError:
                self.cache = {}
            except EOFError:
                self.cache = {}
        else:
            self.cache = {}


    def __call__(self, func):

        @wraps(func)
        def decorated(*args, **kw):
            hasheable = _hash_arg(args)
            if hasheable:
                rtime, result = self.cache.get(hasheable, (None, None))
            else:
                result = None
            if result is not None and time.time() - rtime < self.deadline:
                debug(" Cache load: %s %s %s : %s" % (func.func_name, args,
                    kw, result))
                return result
            else:
                debug(" Cache: No load")
                result = func(*args, **kw)
                if result is not None:
                    hasheable = _hash_arg(args)
                    if hasheable:
                        self.cache[hasheable] = time.time(), result
                        self.count += 1
                        debug(" Cache save: %s %s %s : %s" % (func.func_name,
                            args, kw, result))
                        if self.flush_frequency:
                            if self.count % self.flush_frequency == 0:
                                self.flush()
                else:
                    debug(" %s result is None" % func.func_name)
                return result
        return decorated


    def __del__(self, *args):
        self.flush()


    def flush(self):
        if self.filename:
            file = open(self.filename, "wb")
            pickle.dump(self.cache, file, -1)
            file.close()
            debug("Cache escrito exitosamente en %s" % self.filename)


def main():

    @Cache("fibonar.pickle")
    def fibonar(n):
        if n < 2: return n
        else: return fibonar(n - 1) + fibonar(n - 2)

    print(fibonar(50))
    print(fibonar(250))
    print(fibonar(500))
    print(fibonar(750))
    print(fibonar(1000))

if __name__ == "__main__":
    exit(main())
