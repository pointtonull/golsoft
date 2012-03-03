#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from functools import wraps
import inspect
import os
import pickle
import sys
import time
from debug import error, warning, debug, info
import hashlib

class Zombi:
    def __init__(self):
        self.instances = []
    def append(self, instance):
        self.instances.append(instance)
    def __del__(self):
        for instance in self.instances:
            instance.flush()

_ZOMBI = Zombi()

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
    def __init__(self, filename=None, deadline=0, flush_frequency=0):
        self.count = 0
        self.deadline = deadline
        self.filename = filename
        self.flush_frequency = flush_frequency
        self._ready = False
        self._updated = False
        _ZOMBI.append(self)


    def __call__(self, func):

        @wraps(func)
        def decorated(*args, **kwargs):
            if not self._ready:
                self.load()

            hasheable = _hash_arg(args)
            if hasheable:
                rtime, result = self.cache.get(hasheable, (None, None))
            else:
                result = None

            if result is not None:
                if not self.deadline or time.time() - rtime < self.deadline:
                    debug("Cache load: %s %s %s" % (func.func_name, args,
                        kwargs))
                    return result
                else:
                    debug("Discarting caduced result")
            else:
                debug("Cache: No load")
            result = func(*args, **kwargs)
            if result is not None:
                self.put(args, kwargs, result)
            else:
                debug("%s result is None" % func.func_name)
            return result

        self.func = decorated
        return decorated


    def load(self):
        if self.filename:
            try:
                debug("Opening file cache")
                self.cache = pickle.load(open(self.filename))
            except IOError:
                debug("IOError, creating new empty cache")
                self.cache = {}
            except EOFError:
                debug("EOFError, creating new empty cache2")
                self.cache = {}
        else:
            debug("Creating new empty cache")
            self.cache = {}
        self._ready = True


    def put(self, args, kwargs, result):
        if not self._ready:
            self.load()

        hasheable = _hash_arg(args)
        if hasheable:
            self.cache[hasheable] = time.time(), result
            self.count += 1
            debug("Cache save: %s %s %s" % (self.func.func_name, args, kwargs))
            self._updated = True
            if self.flush_frequency:
                if self.count % self.flush_frequency == 0:
                    debug("Cache autoflushing")
                    self.flush()
        else:
            debug("Cache no hasheable: %s %s %s" % (self.func.func_name,
                args, kwargs))


    def __del__(self, *args):
        self.flush()


    def flush(self):
        if self.filename and self._ready and self._updated:
            file = open(self.filename, "wb")
            pickle.dump(self.cache, file, -1)
            file.close()
            self._updated = False


def main():

    @Cache("fibonar.pickle", 60)
    def fibonar(n):
        if n < 2: return n
        else: return fibonar(n - 1) + fibonar(n - 2)

    print(fibonar(5))
    print(fibonar(50))
    print(fibonar(250))
    print(fibonar(500))
    print(fibonar(750))
    print(fibonar(1000))

if __name__ == "__main__":
    exit(main())
