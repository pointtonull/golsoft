#!/usr/bin/env python
#-*- coding: UTF-8 -*-

import numpy as np

from traits.api import HasTraits, Str, Int, Regex, List, Instance, Enum
from traits.api import on_trait_change
from traits.api import Range
from traitsui.api import View, VGroup, Item, ListEditor

import pea

class Filter(HasTraits):
    name  = Str
    function = Enum(dir(pea),
        label="Function")
    representation = List

    traits_view = View(
        'name',
        'function',
        'representation',
#        width = 0.18,
    )


Cos = Filter(
    name="Cos",
    representation = ["normalize", "equalize", "phase"],
)

pipe = [
    Cos,
]


class ListEditorNotebookSelectionDemo(HasTraits):
    pipe = List(Filter)
    selected = Instance(Filter)
    index = Range(0, 7, mode='spinner')

    traits_view = View(
        Item('index'),
        '_',
        VGroup(
            Item('pipe@',
                  show_label = False,
                  editor = ListEditor(
                      use_notebook = True,
                      deletable = True,
                      selected = 'selected',
                      export = 'DockWindowShell',
                      page_name = '.name')
            )
        ),
    )

    @on_trait_change("selected")
    def update_selector_spin(self, selected):
        self.index = self.pipe.index(selected)

    @on_trait_change("index")
    def update_selected_tab(self, index):
        self.selected = self.pipe[index]


if __name__ == "__main__":
    demo = ListEditorNotebookSelectionDemo(pipe=pipe)
    demo.configure_traits()
