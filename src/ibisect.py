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
        return "Point(%s)" % ", ".join(("%5.3f" % coord
            for coord in self.coordinates))

    def __getattr__(self, attr):
        return getattr(self.coordinates, attr)

    def distance(self, other=None):
        if other is None:
            other = self - self
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
    def __init__(self, function, initial_values=None, plt_args={}):
        dimensions = function.func_code.co_argcount
        if initial_values is None:
            initial_values = []
            for i in range(dimensions + 1):
                point = np.random.normal(0, 1, dimensions)
                distance = np.sum(point ** 2) ** 0.5
                point = (point / distance) * 10
                initial_values.append(point)
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
            plt.axis("off")
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
            print("Left")
            self.returns = -1 
    def key_down(self):
        if self.returns is None:
            print("Draw")
            self.returns = 0
    def key_up(self):
        if self.returns is None:
            print("Draw")
            self.returns = 0
    def key_right(self):
        if self.returns is None:
            print("Right")
            self.returns = 1
    def key_4(self):
        if self.returns is None:
            print("Left")
            self.returns = -1 
    def key_5(self):
        if self.returns is None:
            print("Draw")
            self.returns = 0
    def key_6(self):
        if self.returns is None:
            print("Right")
            self.returns = 1
    def key_q(self):
        if self.returns is None:
            print("Draw")
            self.returns = 0
        plt.close()
        

if __name__ == "__main__":

    from dft import align_phase
    from image import normalize, equalize
    from skimage.filter import canny
    from scipy.ndimage.filters import gaussian_filter
    from autopipe import showimage
    from unwrap import wrapped_gradient, wrapped_diff, unwrap_qg
    import pea

    p = pea.PEA()
    p.filename_holo = "1206-h-det.png"
    p.filename_ref  = "1206-r-det.png"
    p.filename_obj  = "1206-o-det.png"
    p.unvwrapper = unwrap_qg

    def sigmoid(x):
        result = np.abs(1 / (1 + np.exp(-x)))
        return result

#    phase = np.sin(align_phase(p.phase_corrected)[1]) / 2 + 0.5
#    phase = p.unwrapped_phase
#    phase = wrapped_gradient(p.phase_corrected)
#    phase = wrapped_diff(p.phase_corrected)

#    phase = phase / phase.ptp()
#    phase -= phase.min()
    def processor(t_sigma, t_level, sigma, low_threshold, high_threshold):
        t_sigma = sigmoid(t_sigma) * 20
        t_level = sigmoid(t_level)
        sigma = sigmoid(sigma) * 50
        phase = p.unwrapped_phase
        local_context = gaussian_filter(phase, t_sigma)
        phase = (phase - local_context) / (1 - t_level)
        phase /= phase.ptp()
        phase -= phase.min()
        cannied = canny(phase, sigma, low_threshold, high_threshold)
        return np.vstack((phase, cannied))

    initial_values = None
    # Unwrapped
#    initial_values = [(-1.469, -3.429, -4.503, 1.695), (-1.519, -4.142, -4.525, 2.041), (-1.699, -3.178, -4.871, 1.224), (-1.416, -3.636, -5.074, 1.400), (-1.078, -4.102, -4.455, 1.525)]
    # sinoided
#    initial_values = [(2.671, -1.034, -1.752, -11.838), (2.863, -0.547, -1.404, -11.574), (3.073, -0.643, -1.065, -10.938), (3.014, -1.110, -1.634, -11.724), (2.188, -1.027, -1.925, -11.677)]
    # denoised unwrapped
    # [Point(5.095, 4.652, -2.600, -11.167, 0.369), Point(5.050, 4.233, -1.862, -10.801, 0.226), Point(4.725, 4.596, -2.015, -11.064, 0.065), Point(3.842, 4.795, -1.818, -10.624, -0.691), Point(4.410, 4.721, -1.873, -11.075, 0.000), Point(4.824, 4.758, -1.950, -11.196, -0.007)]
    # denoised tonemaped unwrapped phase
    # [Point(11.302, 32.193, -17.454, 0.062, -23.372), Point(11.606, 33.481, -17.760, 0.066, -24.366), Point(10.540, 31.750, -17.362, 0.080, -23.049), Point(9.432, 33.442, -19.050, 0.111, -23.610), Point(11.899, 32.981, -18.453, 0.123, -23.621), Point(11.175, 32.647, -18.048, 0.390, -23.549)]

 
    plt_args = {"cmap":plt.get_cmap("bone")}
    amoeba = Amoeba(processor, initial_values, plt_args=plt_args)
    amoeba.iterate(distance=1)
    print(amoeba.points)
    showimage(amoeba.points[0].value)
