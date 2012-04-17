#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from dft import get_shifted_idft, get_shifted_dft
from image import equalize, imread, imwrite, normalize, get_intensity
from pea import calculate_director_cosines, get_ref_beam, get_auto_mask

from mayavi.core.api import PipelineBase
from mayavi.core.ui.api import SceneEditor
from mayavi.core.ui.mayavi_scene import MayaviScene
from mayavi.tools.mlab_scene_model import MlabSceneModel
from traits.api import HasTraits, Button, File, Range, Enum, Instance, on_trait_change
from traitsui.api import View, Item, Group
from traitsui.menu import OKButton

import scipy
import numpy as np
import sys


class PEA(HasTraits):


    def __init__(self, initialfile=""):
        HasTraits.__init__(self)

        if initialfile:
            self.filename = initialfile
            self.update_image
        self.calculate_ref_beam()



    ## INPUT DATA ##
    filename = File(filter=["*.raw"])
    hologram = scipy.lena()
    grp_datainput = Group(
            Item('filename'),
            label="Input file",
            show_border=True,
        )
    scn_hologram = Instance(MlabSceneModel, ())
    plt_hologram = Instance(PipelineBase)
    vis_inputfile = Item('scn_hologram', editor=SceneEditor(scene_class=MayaviScene),
        height=600, width=600, show_label=False),


    @on_trait_change("filename")
    def update_image(self):
        self.hologram = imread(self.filename)

        if self.plt_hologram is None:
            self.plt_hologram = self.scn_hologram.mlab.imshow(
                self.hologram, colormap="spectral",
                figure=self.scn_hologram.mayavi_scene)
        else:
            self.plt_hologram.mlab_source.set(scalars=self.hologram)



    ## REF-BEAM ##
    ref_beam = np.zeros_like(hologram)
    cos_alpha = Range(-1., 1., 0., mode="xslider")
    cos_beta = Range(-1., 1., 0., mode="xslider")
    btn_director_cosines = Button("Calculate best values")
    ref_beam_vismode = Enum('ref_beam', 'hologram x ref_beam', label="Visualize")
    grp_ref_beam_parameters = Group(
        "cos_alpha",
        "cos_beta",
        Item('btn_director_cosines', show_label=False),
        label="Parameters",
        show_border=True,
    )
    scn_ref_beam = Instance(MlabSceneModel, ())
    plt_ref_beam = Instance(PipelineBase)
    vis_ref_beam = Item('scn_ref_beam', editor=SceneEditor(scene_class=MayaviScene),
        height=600, width=600, show_label=False)


    @on_trait_change("btn_director_cosines")
    def calculate_director_cosines(self):
        cos_alpha, cos_beta = calculate_director_cosines(self.hologram)
        self.cos_alpha = cos_alpha
        self.cos_beta = cos_beta


    @on_trait_change("cos_alpha,cos_beta, ref_beam_vismode")
    def calculate_ref_beam(self):
        self.ref_beam = get_ref_beam(self.hologram.shape, self.cos_alpha, self.cos_beta)
        self.ref_beam_x_hologram = self.ref_beam * self.hologram
        self.spectrum = get_shifted_dft(self.ref_beam_x_hologram)

        if "hologram x"  in self.ref_beam_vismode:
            array = normalize(self.ref_beam_x_hologram)
        else:
            array = self.ref_beam.real

        if self.plt_ref_beam is None:
            self.plt_ref_beam = self.scn_ref_beam.mlab.imshow(
                array, colormap="black-white",
                figure=self.scn_ref_beam.mayavi_scene)
        else:
            self.plt_ref_beam.mlab_source.set(scalars=array)



    ## SPECTRUM MASK ##
    mask = np.zeros_like(hologram)
    softness = Range(0., 30., 0., mode="xslider")
    radious_scale = Range(0., 2., 1., mode="xslider")
    zero_scale = Range(0., 2., 1., mode="xslider")
    cuttop = Range(0., 1., .5, mode="xslider")
    btn_draw_mask = Button("Draw the mask")
    mask_vismode = Enum("hibryd", "mask", "spectrum x mask", label="Visualize")
    grp_mask_parameters = Group(
        "softness",
        "radious_scale",
        "zero_scale",
        "cuttop",
        Item('btn_draw_mask', show_label=False),
        label="Parameters",
        show_border=True,
    )
    scn_mask = Instance(MlabSceneModel, ())
    plt_mask = Instance(PipelineBase)
    vis_mask = Item('scn_mask', editor=SceneEditor(scene_class=MayaviScene),
        height=600, width=600, show_label=False, resizable=True)

    
    @on_trait_change("btn_draw_mask, mask_vismode")
    def generate_mask(self):
        self.spectrum_intensity = get_intensity(self.spectrum)
        self.mask, masked_intensity = get_auto_mask(self.spectrum_intensity,
            self.softness, self.radious_scale, self.zero_scale, self.cuttop)

        if self.mask_vismode == "mask":
            array = self.mask
        elif self.mask_vismode == "hibryd":
            array = normalize(self.mask) + equalize(self.spectrum_intensity)
        else:
            array = equalize(masked_intensity)


        if self.plt_mask is None:
            self.plt_mask = self.scn_mask.mlab.imshow(
                array, colormap="black-white",
                figure=self.scn_mask.mayavi_scene)
        else:
            self.plt_mask.mlab_source.set(scalars=array)



    ## PUT ALL-TOGHETER ##
    view = View(

        Group(
            grp_datainput,
            vis_inputfile,
            label="Data input",
        ),

        Group(
            grp_ref_beam_parameters,
            vis_ref_beam,
            Item('ref_beam_vismode', style='simple', show_label=False),
            label="Reference beam",
        ),

        Group(
            grp_mask_parameters,
            vis_mask,
            Item('mask_vismode', style='simple', show_label=False),
            label="Mask",
        ),

        title="PEA processor",
        buttons=[OKButton], 
        resizable=True,
    )


def main():
    filenames = [filename for filename in sys.argv[1:]]
    filenames = filenames or [""]
    for filename in filenames:
        window = PEA(filename)
        window.configure_traits()

if __name__ == "__main__":
    exit(main())
