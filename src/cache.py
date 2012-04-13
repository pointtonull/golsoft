
#-*- coding: UTF-8 -*-

from debug import error, warning, debug
from functools import wraps
import hashlib
import os
import pickle
import time


HOME = os.path.expanduser("~")
DIRNAME = os.path.join(HOME, ".pycache")
try:
    os.makedirs(DIRNAME)
except OSError, message:
    pass


class ThinDict(dict):
    """
    Similar to builtins dicts except that save hash(key) insteat of key to
    reduce mem/disk consumtion.
    """

    def __setitem__(self, key, value):
        """
        x.__setitem__(key, value) <==> x[key] = value
        """
        key_hash = hash(key)
        debug("ThinDict.__setitem__.key_hash:: %s" % key_hash)
        return dict.__setitem__(self, key_hash, value)


    def __getitem__(self, key):
        """
        x.__getitem__(key) <==> x[key]
        """
        key_hash = hash(key)
        debug("ThinDict.__getitem__.key_hash:: %s" % key_hash)
        return dict.__getitem__(self, key_hash)


    def get(self, key, default):
        """
        D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None.
        """
        key_hash = hash(key)
        if key_hash in self:
            return self.__getitem__(key)
        else:
            return default



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
            with open(self.filename, "wb") as file:
                pickle.dump(object, file, -1)


    def load(self):
        try:
            object = pickle.load(open(self.filename, "rb"))
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
        hash(kwargs)
        return kwargs
    except TypeError:
        try:
            ahash = hashlib.sha1(kwargs).hexdigest()
            return ahash
        except TypeError:
            values = tuple((_hash_args(kwarg) for kwarg in kwargs))
            return values



def make_hasheable(args, kwargs):
    hasheable_args = _hash_args(args)
    hasheable_kwargs = _hash_kwargs(kwargs)
    if hasheable_args is not None and hasheable_kwargs is not None:
        return (hasheable_args, hasheable_kwargs)
    else:
        debug("make_hasheable: %s" % str((hasheable_args, hasheable_kwargs)))
        return None



class Cache:
    def __init__(self, func, ramratio=.5, diskratio=5, deadline=0,
        flush_freq=0, reset=False):
        """
        Diskratio and ramratio are memsize/cputime on MiBs/Secs. A result  
        will be keeped only if size(result) <= cputime * ratio. If a ratio
        value is False this medium will be not used.

        eg:
          todisk = Cache(func, ramratio=False, diskratio=7)
          toram = Cache(func, ramratio=1, diskratio=False)
        """
        debug("Instanced Cache(%s, %s, %s, %s, %s)" % (func, ramratio,
            diskratio, deadline, flush_freq))
        self.func = func
        self.count = 0
        self.deadline = deadline
        self.flush_freq = flush_freq
        self._ready = False
        self._updated = False
        if diskratio is False:
            self.filename = False
        else:
            self.filename = os.path.join(DIRNAME,
                str(hash(func.func_code)) + ".index")
            if reset and os.path.exists(self.filename):
                os.remove(self.filename)
        self.wrapped = self.get_wrapped()
        _ZOMBIE.append(self)


    def __call__(self, *args, **kwargs):
        result = self.wrapped(*args, **kwargs)
        return result


    def get_wrapped(self):

        @wraps(self.func)
        def wrapped(*args, **kwargs):
            if not self._ready:
                self.load()

            hasheable = make_hasheable(args, kwargs)
            if hasheable is not None:
                rtime, pickled = self.cache.get(hasheable, (None, None))
            else:
                pickled = None

            if pickled is not None:
                if not self.deadline or time.time() - rtime < self.deadline:
                    result = pickled.load()
                    debug("Cache load: %s %s -> %s" % (args, kwargs, result))
                else:
                    debug("Discarting caduced result")
                    result = None
            else:
                debug("Cache: No load: %s %s" % (args, kwargs))
                result = None

            if result is None:
                result = self.func(*args, **kwargs)
                if result is not None:
                    self.put(args, kwargs, result)
            return result

        return wrapped


    def load(self):
        if self.filename:
            try:
                debug("Opening index")
                self.cache = pickle.load(open(self.filename, "rb"))
            except IOError:
                debug("IOError, creating new empty index")
                self.cache = ThinDict()
            except EOFError:
                debug("EOFError, creating new empty index")
                self.cache = ThinDict()
        else:
            debug("Creating new empty index")
            self.cache = ThinDict()
        self._ready = True


    def put(self, args, kwargs, result):
        if not self._ready:
            self.load()

        hasheable = make_hasheable(args, kwargs)
        if hasheable is not None:
            self.cache[hasheable] = time.time(), Pickled(result)
            self.count += 1
            debug("Cache save: %s %s %s" % (self.func.func_name, args, kwargs))
            self._updated = True
            if self.flush_freq:
                if self.count % self.flush_freq == 0:
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



class Configurer:
    def __init__(self, **kwargs):
        debug("Instanced Configurer(%s)" % kwargs)
        self.kwargs = kwargs

    def __call__(self, func=None, **kwargs):
        if func:
            wrapped = Cache(func, **self.kwargs).get_wrapped()
            return wrapped
        else:
            self.kwargs.update(kwargs)
            return self


toram = Configurer(diskratio=False)
todisk = Configurer(ramratio=False)
hybrid = Configurer()


def main():

    debug("Main routine")
    @toram
    def fibonar(n):
        if n < 2: return n
        else: return fibonar(n - 1) + fibonar(n - 2)

    print(fibonar(5))
    print(fibonar(250))
    print(fibonar(500))
    print(fibonar(750))

if __name__ == "__main__":
    exit(main())
