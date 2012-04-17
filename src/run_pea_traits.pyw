#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from image import equalize, imread, imwrite, normalize
from mayavi import mlab
from mayavi.core.ui.api import SceneEditor
from mayavi.modules.api import IsoSurface
from mayavi.sources.api import ArraySource
from mayavi.tools.mlab_scene_model import MlabSceneModel
from mayavi.core.ui.mayavi_scene import MayaviScene
from pea import calculate_director_cosines, get_ref_beam
from traits.api import HasTraits, Instance, on_trait_change
from traits.api import HasTraits, Int, Button, File, Range
from traitsui.api import View, Item, Group
from traitsui.menu import OKButton
from mayavi.core.api import PipelineBase
import numpy as np
import traits
import scipy
import sys



from traits.api import HasTraits, Str, Trait
from traitsui.api import Item, Group, View, ImageEnumEditor

image_list = ['ref_beam', 'spectrum x ref_beam']

class Dummy ( HasTraits ):
    """ Dummy class for ImageEnumEditor
    """
    x = Str
    view = View()


class ImageEnumEditorDemo ( HasTraits ):
    """ Defines the ImageEnumEditor demo class.
    """
    image_from_list = Trait(editor=ImageEnumEditor(values = image_list,
        prefix = '@icons:', suffix='_origin', cols=4, klass=Dummy),
        *image_list)
    img_group = Group(
        Item( 'image_from_list', style = 'simple',   label = 'Simple' ),
        Item( '_' ),
        Item( 'image_from_list', style = 'custom',   label = 'Custom' ),
        Item( '_' ),
        Item( 'image_from_list', style = 'text',     label = 'Text' ),
        Item( '_' ),
        Item( 'image_from_list', style = 'readonly', label = 'ReadOnly' )
    )
    view = View(
        img_group,
        title     = 'ImageEnumEditor',
        buttons   = ['OK'],
        resizable = True
    )

ImageEnumEditor().configure_traits()

class PEA(HasTraits):

    filename = File(filter=["*.raw"])
    hologram = scipy.lena()
    scene_hologram = Instance(MlabSceneModel, ())
    plot_hologram = Instance(PipelineBase)

    cos_alpha = Range(-1., 1., 0., mode="xslider")
    cos_beta = Range(-1., 1., 0., mode="xslider")
    btn_director_cosines = Button("Calculate best values")
    scene_ref_beam = Instance(MlabSceneModel, ())
    plot_ref_beam = Instance(PipelineBase)



    def __init__(self, initialfile=""):
        HasTraits.__init__(self)
        self.plot_hologram = self.scene_hologram.mlab.imshow(self.hologram,
            colormap="spectral", figure=self.scene_hologram.mayavi_scene)
        self.plot_ref_beam = self.scene_ref_beam.mlab.imshow(self.hologram,
            colormap="black-white", figure=self.scene_ref_beam.mayavi_scene)
        self.calculate_ref_beam()
        if initialfile:
            self.filename = initialfile
            self.update_image


    @on_trait_change("btn_director_cosines")
    def calculate_director_cosines(self):
        cos_alpha, cos_beta = calculate_director_cosines(self.hologram)
        self.cos_alpha = cos_alpha
        self.cos_beta = cos_beta


    @on_trait_change("filename")
    def update_image(self):
        self.hologram = imread(self.filename)
        self.plot_hologram.mlab_source.set(scalars=self.hologram)


    @on_trait_change("cos_alpha,cos_beta")
    def calculate_ref_beam(self):
        self.ref_beam = get_ref_beam(self.hologram.shape, self.cos_alpha, self.cos_beta)
        self.plot_ref_beam.mlab_source.set(scalars=self.ref_beam.real)


    view = View(
        Group(
            Group(
                Item('filename'),
                label = "Input data",
                show_border = True,
            ),
            Item('scene_hologram', editor=SceneEditor(scene_class=MayaviScene),
                height=600, width=600, show_label=False),
            label = "Data file",
        ),

        Group(
            Group(
                "cos_alpha",
                "cos_beta",
                Item('btn_director_cosines', show_label=False),
                label = "Parameters",
                show_border = True,
            ),
            Item('scene_ref_beam', editor=SceneEditor(scene_class=MayaviScene),
                height=600, width=600, show_label=False),
            label = "Reference beam",
        ),

        title="PEA processor",
        buttons = [OKButton], 
        resizable = True,
    )


def main():
    filenames = [filename for filename in sys.argv[1:]]
    filenames = filenames or [""]
    for filename in filenames:
        window = PEA(filename)
        window.configure_traits()

if __name__ == "__main__":
    exit(main())
