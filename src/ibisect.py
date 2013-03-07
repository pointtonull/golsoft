#!/usr/bin/env python
#-*- coding: UTF-8 -*-

"""
"""

import matplotlib.pyplot as plt
import numpy as np
import cache
from bisect import insort

class Bisector:
    def __init__(self, no, yes=None, min_step=1., search_yes=False):
        self.no = no
        if yes is None:
            self.yes = no + 1
            self.searching_yes = True
        else:
            self.yes = yes
            self.searching_yes = search_yes
        self.min_step = min_step


    def __repr__(self):
        return "%f - %f - %f  %s" % (self.no, self.get_value(), self.yes,
            "(closed)" * self.is_closed())

    def is_closed(self):
        return abs(self.yes - self.no) <= self.min_step

    def get_value(self):
        return (self.no + self.yes) / 2.

    def get_status(self):
        return self.no, self.yes, self.min_step, self.searching_yes

    def set_status(self, status):
        self.no, self.yes, self.min_step, self.searching_yes = status

    def mark_no(self):
        new_no = self.get_value()
        if self.searching_yes:
            self.yes = self.yes * 2 - self.no
        self.no = new_no

    def mark_yes(self):
        if self.searching_yes:
            self.searching_yes = False
        self.yes = self.get_value()


class Point(object):
    def __init__(self, coordinates, function, custom_cmp=cmp):
        self.coordinates = np.array(coordinates)
        self.function = function
        self.cmp = custom_cmp
        self.value = None

    def get_value(self):
        if self.value is None:
            self.value = self.function(*self.coordinates)
        return self.value

    def __cmp__(self, other):
        print("%s vs %s" % (self, other))
        return self.cmp(self.get_value(), other.get_value())

    def __div__(self, other):
        return Point(self.coordinates / other, self.function, self.cmp)

    def __sub__(self, other):
        return Point(self.coordinates - other, self.function, self.cmp)

    def __add__(self, other):
        return Point(self.coordinates + other, self.function, self.cmp)

    def __mul__(self, other):
        return Point(self.coordinates * other, self.function, self.cmp)

    def __repr__(self):
        return "Point(%s)" % ", ".join(("%5.3f" % coord for coord in self.coordinates))

    def __getattr__(self, attr):
        return getattr(self.coordinates, attr)

    def distance(self, other):
        return osum((self - other).coordinates ** 2) ** 0.5



class Handler:
    def __init__(self, updater, bisector, plt_args):
        self.figure = plt.figure()
        self.figure.canvas.mpl_connect('key_press_event', self)
        self.bisector = bisector
        self.updater = updater
        self.plt_args = plt_args
        self.image = plt.imshow(self.updater(self.bisector.get_value()), **plt_args)
        self.memento = [[], []]
        plt.show()

    def __call__(self, event):
        if not event.key is None:
            attr = getattr(self, "key_" + event.key, None)
            if attr:
                attr()

    def key_down(self):
        """Mark the current value as no"""
        self.bisector.mark_no()
        self.update_status()

    def update_status(self):
        if self.bisector.is_closed():
            self.quit()
            print(self.bisector)
        else:
            self.memento[0].append(self.bisector.get_status())
            self.memento[1] = []
            print(self.bisector)
            self.update_image()

    def key_up(self):
        """Mark the current value as yes"""
        self.bisector.mark_yes()
        self.update_status()

    def key_left(self):
        """Rudimentary undo feature"""
        if self.memento[0]:
            self.memento[1].append(self.bisector.get_status())
            self.bisector.set_status(self.memento[0].pop())
            print(self.bisector)
            self.update_image()

    def key_right(self):
        """Rudimentary redo feature"""
        if self.memento[1]:
            self.memento[0].append(self.bisector.get_status())
            self.bisector.set_status(self.memento[1].pop())
            print(self.bisector)
            self.update_image()

    def key_q(self):
        self.quit()

    def quit(self):
        """Close the search now"""
        print("Quit")
        plt.close()

    def update_image(self):
        """Render the image with the new value"""
        self.image.set_array(self.updater(self.bisector.get_value()))
        self.figure.canvas.draw()
        print("Image updated")


class Comparable(object):
    def __init__(self, thing, cmp_fnc):
        self.thing = thing
        self.cmp = cmp_fnc

    def __cmp__(self, other_thing):
        return self.cmp(self, other_thing)

    def __getattr__(self, attr):
        return getattr(self.thing, attr)


def osum(iterable):
    initial = 0
    for item in iterable:
        if initial is 0:
            initial = item
        else:
            initial = initial + item
    return initial

def omean(iterable):
    initial = 0
    lenght = 0
    for item in iterable:
        lenght += 1
        if initial is 0:
            initial = item
        else:
            initial = initial + item
    return initial / float(lenght)

class Amoeba:
    def __init__(self, function, initial_values=None, plt_args=None):
        dimensions = function.func_code.co_argcount
        if initial_values is None:
            initial_values = [np.random.normal(0, 15, dimensions)
                for i in range(dimensions + 1)]            
        else:
            assert len(initial_values) == dimensions + 1

        self.ndim = dimensions
        self.cmp = Icmp(plt_args)
        self.points = []
        for coordinates in initial_values:
            insort(self.points, Point(coordinates, function, self.cmp))


    def iterate(self, iterations=None, distance=None):
        if iterations:
            for iteraction in range(iterations):
                self.iterate()
        elif distance:
            while self.points[0].distance(self.points[-1]) > distance:
                self.iterate()
        else:
            # reflection
            assert len(self.points) == self.ndim + 1
            x_o = omean(self.points[:-1])
            x_r = x_o * 2 - self.points[-1]
            if x_r < self.points[0]:
                print("expansion")
                x_e = x_r * 2 - x_o
                if x_e < x_r:
                    self.points.pop()
                    self.points.insert(0, x_e)
                    return
                else:
                    self.points.pop()
                    self.points.insert(0, x_r)
                    return
            elif x_r < self.points[-2]:
                print("reflection")
                self.points.pop()
                insort(self.points, x_r)
                return
            else:
                x_c = (x_o + self.points[-1]) / 2.
                if x_c < self.points[-1]:
                    print("contraction")
                    self.points.pop()
                    insort(self.points, x_c)
                else:
                    print("reduction")
                    newpoints = []
                    for point in self.points:
                        insort(newpoints, (point + self.points[0]) / 2.)
                    self.points = newpoints
        

class Icmp:
    def __init__(self, plt_args={}):
        self.figure = plt.figure()
        self.figure.canvas.mpl_connect('key_press_event', self.handler)
        self.plt_args = plt_args
        self.image = None
        self.returns = None

    @cache.toram
    def __call__(self, left, right):
        if (left == right).all():
            print("Auto draw!")
            return 0
        image = np.hstack((left, right))
        if self.image is None:
            self.image = plt.imshow((image), **self.plt_args)
            plt.show(block=False)
        else:
            self.image.set_array(image)
            self.image.set_visible(True)
            self.figure.canvas.draw()

        while self.returns is None:
            plt.waitforbuttonpress()
        result = self.returns
        self.image.set_visible(False)
        self.figure.canvas.draw()
        self.returns = None
        return result

    def handler(self, event):
        if not event.key is None:
            attr = getattr(self, "key_" + event.key, None)
            if attr:
                attr()

    def key_left(self):
        if self.returns is None:
            self.returns = -1 
    def key_down(self):
        if self.returns is None:
            self.returns = 0
    def key_up(self):
        if self.returns is None:
            self.returns = 0
    def key_right(self):
        if self.returns is None:
            self.returns = 1
    def key_4(self):
        if self.returns is None:
            self.returns = -1 
    def key_5(self):
        if self.returns is None:
            self.returns = 0
    def key_6(self):
        if self.returns is None:
            self.returns = 1
    def key_q(self):
        if self.returns is None:
            self.returns = 0
        plt.close()
        

if __name__ == "__main__":

#    from dft import align_phase
#    import pea
#    import unwrap

#    tau = np.pi * 2

#    p = pea.PEA()
#    p.unwrapper = unwrap.unwrap_qg
#    p.filename_holo = "1206-h-det.png"
#    p.filename_ref  = "1206-r-det.png"
#    p.filename_obj  = "1206-o-det.png"

#    def updater(distance):
#        print("Distance: %f" % distance)
#        p.distance = distance
#        return align_phase(p.phase)[1]
            
#    plt_args = {"cmap":plt.get_cmap("bone")}
#    handler = Handler(updater, Bisector(-5, 5, 0.0001, search_yes=True), plt_args=plt_args)
#    p.phase_denoise = handler.bisector.yes


    from scipy.misc import lena
    from image import normalize, equalize
    from autopipe import showimage
    from scipy.ndimage.filters import median_filter, gaussian_filter
    from dft import align_phase

#    icmp = Icmp({"cmap":plt.get_cmap("bone")})
#    lenas = sortedlist()
#    for wrapp in np.random.random_integers(0, 255, 10):
#        lenas.add(Comparable(align_phase(normalize(lena() % wrapp) / 40.74)[1], icmp))
#    showimage(np.hstack(lenas))

    def sigmoid(x):
        result = np.abs(1 / (1 + np.exp(-x)))
        return result

    image = lena()
    def processor(median_size, tonemapping_level, equalize_level):
        median_size = sigmoid(median_size) * 50
        tonemapping_level = sigmoid(tonemapping_level)
        equalize_level = sigmoid(equalize_level)
        local_context = gaussian_filter(image, median_size)
        tonemapped = normalize(image - local_context * tonemapping_level)
        equalized = tonemapped * (1 - equalize_level) + equalize(tonemapped) * equalize_level
        return equalized
    plt_args = {"cmap":plt.get_cmap("bone")}
    amoeba = Amoeba(processor, plt_args=plt_args)
    amoeba.iterate(distance=1)
    print(amoeba.points)
    showimage(np.hstack((image, amoeba.points[0].value)))
