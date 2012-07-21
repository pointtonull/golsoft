#!/usr/bin/env python
#-*- coding: UTF-8 -*-

import sys

import numpy as np
import scipy

from mayavi.core.api import PipelineBase
from mayavi.core.ui.api import SceneEditor
from mayavi.core.ui.mayavi_scene import MayaviScene
from mayavi.tools.mlab_scene_model import MlabSceneModel
from traits.api import HasTraits, Button, File, Range, Enum, Instance, Bool
from traits.api import on_trait_change
from traitsui.api import View, Item, Group, HSplit
from traitsui.menu import OKButton

from dft import get_shifted_idft, get_shifted_dft
from image import equalize, imread, imwrite, normalize, get_intensity
from pea import calculate_director_cosines, get_ref_beam, get_auto_mask
from pea import get_propagation_array


class PEA(HasTraits):

    def __init__(self, initialfile=""):
        HasTraits.__init__(self)

        if initialfile:
            self.filename = initialfile
            self.update_image
        self.calculate_ref_beam()
    

    ## OVERVIEW ##

    filename = File(filter=[u"*.raw"])
    overview_vismode = Enum("input map", "input surface", "phase map",
        "phase surface", "module map", "hibryd surface", label="Visualize")

    scn_overview = Instance(MlabSceneModel, ())
    plt_overview = Instance(PipelineBase)
    vis_overview = Item('scn_overview', 
        editor=SceneEditor(scene_class=MayaviScene), height=600, width=600,
        show_label=False),

    grp_datainput = Group(
        Item('filename'),
        label="Input file",
        show_border=True,
    )

    grp_overview_visualizer = Group(
        vis_overview,
        Item('overview_vismode', style='simple', show_label=False),
    )

    @on_trait_change("filename")
    def update_image(self):
        self.hologram = imread(self.filename)

        if self.plt_overview is None:
            self.plt_hologram = self.scn_overview.mlab.imshow(
                self.hologram, colormap="spectral",
                figure=self.scn_overview.mayavi_scene)
        else:
            self.plt_hologram.mlab_source.set(scalars=self.hologram)



    ## REF-BEAM ##
    use_ref_beam = Bool(False)


    cos_alpha = Range(-1., 1., 0., mode="xslider")
    cos_beta = Range(-1., 1., 0., mode="xslider")
    btn_director_cosines = Button("Calculate best values")
    ref_beam_vismode = Enum('ref_beam', 'hologram x ref_beam', 
        label="Visualize")

    scn_ref_beam = Instance(MlabSceneModel, ())
    plt_ref_beam = Instance(PipelineBase)
    vis_ref_beam = Item('scn_ref_beam', 
        editor=SceneEditor(scene_class=MayaviScene), height=600, width=600,
        show_label=False)

    grp_ref_beam_parameters = Group(
        "cos_alpha",
        "cos_beta",
        Item('btn_director_cosines', show_label=False),
        label="Reference beam parameters",
        show_border=True,
        enabled_when="use_ref_beam",
    )

    grp_ref_beam_visualizer = Group(
        vis_ref_beam,
        Item('ref_beam_vismode', style='simple', show_label=False),
    )


    @on_trait_change("btn_director_cosines")
    def calculate_director_cosines(self):
        cos_alpha, cos_beta = calculate_director_cosines(self.hologram)
        self.cos_alpha = cos_alpha
        self.cos_beta = cos_beta


    @on_trait_change("cos_alpha,cos_beta, ref_beam_vismode")
    def calculate_ref_beam(self):
        self.ref_beam = get_ref_beam(self.hologram.shape, self.cos_alpha,
            self.cos_beta)
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
    use_masking = Bool(True)

    softness = Range(0., 30., 0., mode="xslider", enter_set=True,
        auto_set=False)
    radious_scale = Range(0., 2., 1., mode="xslider", enter_set=True,
        auto_set=False)
    use_zero_mask = Bool(True)
    zero_scale = Range(0., 2., 1., mode="xslider", enter_set=True,
        auto_set=False, enabled_when="use_zero_mask")
    use_cuttop = Bool(False)
    cuttop = Range(.99, 1., 0., mode="xslider", enter_set=True,
        auto_set=False, enabled_when="use_cuttop")
    mask_vismode = Enum("hibryd", "mask", "spectrum x mask",
        label="Visualize")

    scn_mask = Instance(MlabSceneModel, ())
    plt_mask = Instance(PipelineBase)
    vis_mask = Item('scn_mask', editor=SceneEditor(scene_class=MayaviScene),
        height=600, width=600, show_label=False, resizable=True)

    grp_mask_parameters = Group(
        "softness",
        "radious_scale",
        "use_zero_mask",
        "zero_scale",
        "use_cuttop",
        "cuttop",
        label="Spectrum mask parameters",
        show_border=True,
        enabled_when="use_masking",
    )

    grp_mask_visualizer = Group(
        vis_mask,
        Item('mask_vismode', style='simple', show_label=False),
    )

    
    @on_trait_change("mask_vismode, softness, radious_scale, zero_scale")
    def generate_mask(self):
        self.spectrum_intensity = get_intensity(self.spectrum)
        self.mask, self.masked_hologram = get_auto_mask(
            self.hologram, self.softness, self.radious_scale,
            self.zero_scale)

        if self.mask_vismode == "mask":
            array = self.mask
        elif self.mask_vismode == "hibryd":
            array = normalize(self.mask) + equalize(self.spectrum_intensity)
        else:
            array = equalize(self.masked_hologram)

        if self.plt_mask is None:
            self.plt_mask = self.scn_mask.mlab.imshow(
                array, colormap="black-white",
                figure=self.scn_mask.mayavi_scene)
        else:
            self.plt_mask.mlab_source.set(scalars=array)



    ## PROPAGATION ##
    use_propagation = Bool(True)

    distance = Range(-15., 15., 0., mode="xslider", enter_set=True,
        auto_set=False)
    propagation_vismode = Enum("module masked", "phase masked",
        "module all", "phase all", label="Visualize")

    scn_propagation = Instance(MlabSceneModel, ())
    plt_propagation = Instance(PipelineBase)
    vis_propagation = Item('scn_propagation',
        editor=SceneEditor(scene_class=MayaviScene), height=600, width=600,
        show_label=False, resizable=True)

    grp_propagation_parameters = Group(
        "distance",
        label="Propagation parameters",
        show_border=True,
        enabled_when="use_propagation",
    )

    grp_propagation_visualizer = Group(
        vis_propagation,
        Item('propagation_vismode', style='simple', show_label=False),
    )

    @on_trait_change("propagation_vismode, distance")
    def propagate_hologram(self):
        self.propagation_array = get_propagation_array(self.hologram.shape,
            self.distance)
        self.propagated = self.propagation_array * self.masked_hologram
    

    ## UNWRAPPING ##
    use_unwrapping = Bool(True)


    def __init__(self, initialfile=""):
        HasTraits.__init__(self)

        if initialfile:
            self.filename = initialfile
            self.update_image
        self.calculate_ref_beam()



    ## PUT ALL-TOGHETER ##
    view = View(

        Group(
            HSplit(
                Group(
                    grp_datainput,
                    "use_ref_beam",
                    grp_ref_beam_parameters,
                    "use_masking",
                    grp_mask_parameters,
                    "use_propagation",
                    grp_propagation_parameters,
                ),
                grp_overview_visualizer,
            ),
            label="Overview",
        ),

        Group(
            HSplit(
                grp_ref_beam_parameters,
                grp_ref_beam_visualizer,
            ),
            label="Reference beam",
        ),

        Group(
            HSplit(
                grp_mask_parameters,
                grp_mask_visualizer,
            ),
            label="Mask",
        ),

        Group(
            HSplit(
                grp_propagation_parameters,
                grp_propagation_visualizer,
            ),
            label="Propagation",
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
