#!/usr/bin/env python
#-*- coding: UTF-8 -*-


from traits.api import HasTraits, Int, Button, Float, File, Range
from traitsui.api import View, Item, Group, ImageEditor
from traitsui.menu import OKButton
from pea import get_pea, calculate_director_cosines
from image import equalize, imread
import Image as pil

class Counter(HasTraits):
    alpha = Range(-1., 1., 0., mode="xslider")
    beta = Range(-1., 1., 0., mode="xslider")
    calculate_angles = Button()
    filename = File(filter=["*.raw"])
    rawimage = None
    bitmap = None
    image = ImageEditor()

    def _calculate_angles_fired(self):
        cos_alpha, cos_beta = calculate_director_cosines(self.rawimage)
        self.alpha = cos_alpha
        self.beta = cos_beta


    def _filename_changed(self):
        print(self.filename)
        self.rawimage = imread(self.filename)
        self.bitmap = pil.fromarray(self.rawimage)
        self.image = ImageEditor(image=self.bitmap)


    view = View(
        Group(
            Group(
                Item('filename', emphasized=True),
                "alpha",
                "beta",
                Item('calculate_angles', show_label=False),
                "image",
                label = "Reference beam",
                show_border = True,
            ),
            label = "Data file",
        ),
        Group(
            Group(
                Item('filename', emphasized=True),
                "alpha",
                "beta",
                Item('calculate_angles', show_label=False),
                "image",
                label = "Reference beam",
                show_border = True,
            ),
            label = "Data file",
        ),
        buttons = [OKButton], 
        resizable = True,
    )


def main():
    counter = Counter()
    counter.configure_traits()

if __name__ == "__main__":
    exit(main())
