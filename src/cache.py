#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from debug import error, warning, debug, info
from functools import wraps
import hashlib
import inspect
import os
import pickle
import sys
import time


HOME = os.path.expanduser("~")
DIRNAME = os.path.join(HOME, ".pycache")
try:
    os.makedirs(DIRNAME)
except OSError, message:
    pass


class Zombie:
    """
    helper to force flush on deleting (instance as a global _* variable)
    """

    def __init__(self):
        self.instances = []

    def append(self, instance):
        self.instances.append(instance)

    def __del__(self):
        for instance in self.instances:
            instance.flush()


_ZOMBIE = Zombie()


class Pickled:
    """
    based on the assumption of that all instances are inmutables respect to 
    their hashes. Is quite reasonable, isnt it?
    """

    def __init__(self, object):
        self.hash = hash(_hash_args(object))
        self.filename = os.path.join(DIRNAME, str(self.hash) + ".pickle")
        if not os.path.exists(self.filename):
            with open(self.filename, "w") as file:
                pickle.dump(object, file, -1)


    def load(self):
        try:
            object = pickle.load(open(self.filename))
        except IOError:
            object = None
        return object



def _hash_args(args):
    try:
        hash(args)
        return args
    except TypeError:
        try:
            ahash = hashlib.sha1(args).hexdigest()
            return ahash
        except TypeError:
            values = tuple((_hash_args(arg) for arg in args))
            return values



def _hash_kwargs(kwargs):
    try:
        hash(args)
        return args
    except TypeError:
        try:
            ahash = hashlib.sha1(args).hexdigest()
            return ahash
        except TypeError:
            values = tuple((_hash_args(arg) for arg in args))
            return values



class Cache:
    def __init__(self, filename=None, deadline=0, flush_frequency=0):
        #TODO: add ratio time/size bound
        self.count = 0
        self.deadline = deadline
        self.filename = filename
        self.flush_frequency = flush_frequency
        self._ready = False
        self._updated = False
        _ZOMBIE.append(self)


    def __call__(self, func):

        @wraps(func)
        def decorated(*args, **kwargs):
            if not self._ready:
                self.load()

            hasheable = _hash_args(args)
            if hasheable:
                rtime, pickled = self.cache.get(hasheable, (None, None))
            else:
                pickled = None

            if pickled is not None:
                if not self.deadline or time.time() - rtime < self.deadline:
                    debug("Cache load: %s %s %s" % (func.func_name, args,
                        kwargs))
                    return pickled.load()
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
        #TODO: must include hash(func.func_code.co_code) as first key value
        #FIXME: must include kawargs as part of key value
        if not self._ready:
            self.load()

        hasheable = _hash_args(args)
        if hasheable:
            self.cache[hasheable] = time.time(), Pickled(result)
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
