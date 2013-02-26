#!/usr/bin/env python
#-*- coding: UTF-8 -*-

"""
"""

import matplotlib.pyplot as plt
import numpy as np


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



class IBisector:
    def __init__(self, updater, bisector, plt_args):
        self.figure = plt.figure()
        self.figure.canvas.mpl_connect('key_press_event', self)
        self.bisector = bisector
        self.updater = updater
        self.plt_args = plt_args
        self.image = plt.imshow(self.updater(self.bisector.get_value()), **plt_args)
        self.memento = [[], []]
        print("Image updated")
        plt.show()

    def __call__(self, event):
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


"""
threshold = 0.01
while f(poins[-1]) - f(points[0] > threshold:

#1. settle
    Ordenar de acuerdo a f(x) los vertices tal que:
    f(points[0]) ≤ f(points[1]) ≤ f(points[2]) ≤ … ≤ f(points[-1])
#2. centroid
    x_o = mean(points[:-1])
#3. reflection
    x_r = 2 * x_o - points[-1]
    if  f(points[0]) ≤ f(x_r) < f(points[-2]):
        points[-1] = x_r
        continue
#4. expansion
    if f(x_r) < f(points[0]):
        x_e = x_o + 2 * (x_o - points[-1])
        if f(x_e) < f(x_r):
            points[-1] = x_e
        else:
            points[-1] = r_r
        continue
#5. contraction
    assert f(x_r) ≥ f(points[-2])
    x_c = x_o - (x_o - points[-1]) / 2
    if f(x_c) < f(points[-1]):
        points[-1] = x_c
        continue
#6. reduction
    for i in range(1, n+1):
        points[i] = points[0] + (points[i] - points[0]) / 2
        continue
"""


if __name__ == "__main__":

    from dft import align_phase
    import pea
    import unwrap

    tau = np.pi * 2

    p = pea.PEA()
    p.unwrapper = unwrap.unwrap_qg
    p.filename_holo = "1206-h-det.png"
    p.filename_ref  = "1206-r-det.png"
    p.filename_obj  = "1206-o-det.png"

    def updater(distance):
        print("Distance: %f" % distance)
        p.distance = distance
        return align_phase(p.phase)[1]
            
    plt_args = {"cmap":plt.get_cmap("bone")}
    handler = Handler(updater, Bisector(-5, 5, 0.0001, search_yes=True), plt_args=plt_args)
    p.phase_denoise = handler.bisector.yes
