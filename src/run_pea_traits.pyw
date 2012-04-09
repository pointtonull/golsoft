#!/usr/bin/env python
#-*- coding: UTF-8 -*-


from traits.api import HasTraits, Int, Button, Float, File, Range
from traitsui.api import View, Item, Group
from traitsui.menu import OKButton
from math import cos
from pea import get_pea, guess_director_angles
from image import equalize, imread, normalize
import random

class Counter(HasTraits):
    cos_alpha =  Float()
    cos_beta = Float()
    calculate_angles = Button()
    filename = File(filter=["*.raw"])
    image = None

    def _calculate_angles_fired(self):
        alpha, beta = guess_director_angles(self.image)
        self.cos_alpha = cos(alpha)
        self.cos_beta = cos(beta)


    def _filename_changed(self):
        print(self.filename)
        self.image = imread(self.filename) 


    view = View(
        Item('filename', emphasized=True),
        Group(
            "alpha",
            "beta",
            Item('calculate_angles', show_label=False),
            label = "Reference beam",
#            show_border = True,
        ),
        buttons = [OKButton], 
        resizable = True,
    )


def main():
    counter = Counter()
    counter.configure_traits()

if __name__ == "__main__":
    exit(main())
