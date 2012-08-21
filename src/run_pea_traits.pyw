#!/usr/bin/env python
#-*- coding: UTF-8 -*-

import sys

from mayavi.core.api import PipelineBase
from mayavi.core.ui.api import SceneEditor
from mayavi.core.ui.mayavi_scene import MayaviScene
from mayavi.tools.mlab_scene_model import MlabSceneModel
from traits.api import Bool, Str
from traits.api import HasTraits, Button, File, Range, Enum, Instance
from traits.api import on_trait_change
from traitsui.api import View, Item, Group, HSplit
from traitsui.menu import OKButton
import numpy as np

from dft import get_shifted_idft, get_shifted_dft
from image import equalize, imread, normalize
from pea import calculate_director_cosines, get_ref_beam, get_auto_mask
from pea import get_propagation_array, get_phase, get_module
from unwrap import unwrap_qg, unwrap_wls
from autofocus import guess_focus_distance


class PEA(HasTraits):

    def __init__(self, initialfile=""):
        HasTraits.__init__(self)

        if initialfile:
            self.filename = initialfile
    


    ## STAGES ##
    empty = np.zeros((200, 200))
    hologram = empty # ready
    ref_beam = empty # ready
    r_hologram = empty # ready
    spectrum = empty # ready
    centered_spectrum = empty # ready
    mask = empty # ready
    masked_spectrum = empty # ready
    propagated = empty # ready
    reconstructed = empty # ready
    module = empty # ready
    wrapped_phase = empty # ready
    unwrapped_phase = empty # ready



    ## OVERVIEW ##
    filename = File(filter=[u"*.raw"])
    overview_vismode = Enum("input map", "spectrum", "module", "phase map",
        "unwrapped phase", "unwrapped phase map", label="Visualize")

    scn_overview = Instance(MlabSceneModel, ())
    plt_overview = None
    plt_overview_surf = None
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
    def update_hologram(self):
        print("Updating hologram")
        self.hologram = imread(self.filename)
        self.empty = np.zeros_like(self.hologram)
        self.update_ref_beam()
        self.update_overview_vis()

    
    @on_trait_change("overview_vismode")
    def update_overview_vis(self):
        vismode = self.overview_vismode
        rep_type = "image"
        if vismode == "input map":
            array = self.hologram
            color = "spectral"
        elif vismode == "spectrum":
            array = normalize(self.mask) + equalize(self.centered_spectrum)
            color = "gist_stern"
        elif vismode == "module":
            array = self.module
            color = "bone"
        elif vismode == "phase map":
            array = self.wrapped_phase
            color = "spectral"
        elif vismode == "unwrapped phase":
            array = self.unwrapped_phase
            color = "bone"
            rep_type = "surface"
        elif vismode == "unwrapped phase map":
            array = self.unwrapped_phase
            color = "spectral"
            
        else:
            print("Unrecognized option on overview_vismode: %s" % 
                vismode)
            array = self.empty
            color = "spectral"


        if rep_type == "image":
            if self.plt_overview_surf:
                self.plt_overview_surf.visible = False
            if self.plt_overview is None:
                print("Creating new graph")

                self.plt_overview = self.scn_overview.mlab.imshow(
                    array, colormap=color,
                    figure=self.scn_overview.mayavi_scene)
            else:
                self.plt_overview.visible = True
                self.plt_overview.mlab_source.lut_type = color
                self.plt_overview.mlab_source.scalars = array
        else:
            warp_scale = -100 / array.ptp()
            if self.plt_overview:
                self.plt_overview.visible = False
            if self.plt_overview_surf is None:
                print("Creating new surface graph")

                self.plt_overview_surf = self.scn_overview.mlab.surf(
                    array, colormap=color, warp_scale=warp_scale,
                    figure=self.scn_overview.mayavi_scene)
            else:
                self.plt_overview_surf.visible = True
                self.plt_overview_surf.mlab_source.lut_type = color
                self.plt_overview_surf.mlab_source.scalars = array
                self.plt_overview_surf.mlab_source.warp_scale = warp_scale



    ## REF-BEAM ##
    use_ref_beam = Bool(False)

    use_auto_angles = Bool(True)
    cos_alpha = Range(-1., 1., 0., mode="xslider", enter_set=True,
        auto_set=False)
    cos_beta = Range(-1., 1., 0., mode="xslider", enter_set=True,
        auto_set=False)
    ref_beam_vismode = Enum('ref_beam', 'hologram x ref_beam', 
        label="Visualize")

    scn_ref_beam = Instance(MlabSceneModel, ())
    plt_ref_beam = Instance(PipelineBase)
    vis_ref_beam = Item('scn_ref_beam', 
        editor=SceneEditor(scene_class=MayaviScene), height=600, width=600,
        show_label=False)

    grp_ref_beam_parameters = Group(
        "use_ref_beam",
        Group(
            "use_auto_angles",
            Group(
                "cos_alpha",
                "cos_beta",
                enabled_when="not use_auto_angles"
            ),
            label="Reference beam parameters",
            show_border=True,
            enabled_when="use_ref_beam",
        )
    )

    grp_ref_beam_visualizer = Group(
        vis_ref_beam,
        Item('ref_beam_vismode', style='simple', show_label=False),
    )


    @on_trait_change("use_ref_beam, use_auto_angles, hologram")
    def calculate_director_cosines(self):
        if self.use_ref_beam and self.use_auto_angles:
            cos_alpha, cos_beta = calculate_director_cosines(self.hologram)
            self.cos_alpha, self.cos_beta = cos_alpha, cos_beta


    @on_trait_change("cos_alpha, cos_beta, use_ref_beam")
    def update_ref_beam(self):
        if self.use_ref_beam:
            print("Updating reference beam")
            self.ref_beam = get_ref_beam(self.hologram.shape, self.cos_alpha,
                self.cos_beta)
            self.r_hologram = self.ref_beam * self.hologram
        else:
            print("Not updating reference beam")
            self.ref_beam = self.empty
            self.r_hologram = self.hologram

        self.spectrum = get_shifted_dft(self.r_hologram)

        self.update_ref_beam_vis()
        self.update_overview_vis()
        self.update_mask()


    @on_trait_change("ref_beam_vismode")
    def update_ref_beam_vis(self):
        if "hologram x"  in self.ref_beam_vismode:
            array = get_module(self.r_hologram)
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
    zero_scale = Range(0., 2., 1.4, mode="xslider", enter_set=True,
        auto_set=False, enabled_when="use_zero_mask")
    use_cuttop = Bool(False)
    cuttop = Range(.00, .01, 0., mode="xslider", enter_set=True,
        auto_set=False, enabled_when="use_cuttop")
    mask_vismode = Enum("hibryd", "mask", "spectrum",
        label="Visualize")

    scn_mask = Instance(MlabSceneModel, ())
    plt_mask = Instance(PipelineBase)
    vis_mask = Item('scn_mask', editor=SceneEditor(scene_class=MayaviScene),
        height=600, width=600, show_label=False, resizable=True)

    grp_mask_parameters = Group(
        "use_masking",
        Group(
            "softness",
            "radious_scale",
            "use_zero_mask",
            Item("zero_scale", enabled_when="use_zero_mask"),
            "use_cuttop",
            Item("cuttop", enabled_when="use_cuttop"),
            label="Spectrum mask parameters",
            show_border=True,
            enabled_when="use_masking",
        )
    )

    grp_mask_visualizer = Group(
        vis_mask,
        Item('mask_vismode', style='simple', show_label=False),
    )

    
    @on_trait_change("use_masking, softness, radious_scale, use_zero_mask,"
        "zero_scale, use_cuttop, cuttop")
    def update_mask(self):
        print("Updating mask")
        if self.use_masking:
            zero_scale = self.zero_scale if self.use_zero_mask else 0
            cuttop = self.cuttop if self.use_cuttop else 0
            mask, masked_spectrum, centered_spectrum = get_auto_mask(
                self.spectrum, self.softness, self.radious_scale,
                zero_scale, cuttop)
            self.mask = mask
            self.masked_spectrum = masked_spectrum
            self.centered_spectrum = centered_spectrum

            self.update_mask_vis()
            self.update_overview_vis()
        else:
            print("Not using mask")

        self.update_propagation()


    def update_mask_vis(self):
        if self.mask_vismode == "hibryd":
            array = normalize(self.mask) + equalize(self.centered_spectrum)
        elif self.mask_vismode == "mask":
            array = normalize(equalize(self.masked_spectrum))
        elif self.mask_vismode == "spectrum":
            array = equalize(self.centered_spectrum)
        else:
            print("Unrecognized option on mask_vismode: %s" % 
                self.mask_vismode)

        if self.plt_mask is None:
            self.plt_mask = self.scn_mask.mlab.imshow(
                array, colormap="gist_stern",
                figure=self.scn_mask.mayavi_scene)
        else:
            self.plt_mask.mlab_source.set(scalars=array)


    ## PROPAGATION ##
    use_propagation = Bool(True)

    btn_guess_focus = Button("Guess focus distance")
    distance = Range(-0.20, 0.20, 0., mode="xslider", enter_set=True,
        auto_set=False)
    propagation_vismode = Enum("module", "phase", label="Visualize")

    scn_propagation = Instance(MlabSceneModel, ())
    plt_propagation = Instance(PipelineBase)
    vis_propagation = Item('scn_propagation',
        editor=SceneEditor(scene_class=MayaviScene), height=600, width=600,
        show_label=False, resizable=True)

    grp_propagation_parameters = Group(
        "use_propagation",
        Group(
            "btn_guess_focus",
            "distance",
            label="Propagation parameters",
            show_border=True,
            enabled_when="use_propagation",
        )
    )

    grp_propagation_visualizer = Group(
        vis_propagation,
        Item('propagation_vismode', style='simple', show_label=False),
    )


    @on_trait_change("btn_guess_focus")
    def guess_focus_distance(self):
        self.distance = guess_focus_distance(self.masked_spectrum)
    

    @on_trait_change("propagation_vismode, distance")
    def update_propagation(self):
        propagation_array = get_propagation_array(self.hologram.shape,
            self.distance)
        self.propagated = propagation_array * self.masked_spectrum
        self.reconstructed = get_shifted_idft(self.propagated)
        self.module = get_module(self.reconstructed)
        self.wrapped_phase = get_phase(self.reconstructed)
        self.update_propagation_vis()
        self.update_overview_vis()
        self.update_unwrapping_phase()


    @on_trait_change("propagation_vismode")
    def update_propagation_vis(self):
        if self.propagation_vismode == "module":
            array = self.module
            color = "bone"
        elif self.propagation_vismode == "phase":
            array = self.wrapped_phase
            color = "bone"
        else:
            print("Unrecognized option on propagation_vismode: %s" % 
                self.mask_vismode)
            color = "prism"

        if self.plt_propagation is None:
            self.plt_propagation = self.scn_propagation.mlab.imshow(
                array, colormap=color,
                figure=self.scn_propagation.mayavi_scene)
        else:
            self.plt_propagation.mlab_source.set(scalars=array)


    ## UNWRAPPING ##
    use_unwrapping = Bool(True)

    unwrapping_method = Enum("Unweighted Least Squares", "Quality Guided")
    unwrapping_vismode = Enum("phase", "hibryd", label="Visualize")
    txt_unwrapping_method = Str(
        """
Unweighted Least Squares: 
Uses a simple wavelet aproach on with applys then apply a least squares
regression to create a fast but inexact approximation.

Quality Guided:
Uses the module image as a quality map to determine the best
path to apply the unwrapping""")

    scn_unwrapping = Instance(MlabSceneModel, ())
    plt_unwrapping = Instance(PipelineBase)
    vis_unwrapping = Item('scn_unwrapping',
        editor=SceneEditor(scene_class=MayaviScene), height=600, width=600,
        show_label=False, resizable=True)

    grp_unwrapping_parameters = Group(
        "use_unwrapping",
        Group(
            Item("unwrapping_method", show_label=False),
            Item("txt_unwrapping_method", style='readonly', show_label=False),
            label="Unwrapping method",
            show_border=True,
            enabled_when="use_unwrapping",
        )
    )

    grp_unwrapping_visualizer = Group(
        vis_unwrapping,
        Item('unwrapping_vismode', style='simple', show_label=False),
    )


    @on_trait_change("unwrapping_method")
    def update_unwrapping_phase(self):
        if self.use_unwrapping:
            print("Unwrapping phase")
            method = self.unwrapping_method
            if "Unweighted Least Squares" in method:
                self.unwrapped_phase = unwrap_wls(self.wrapped_phase)
            elif "Quality Guided" in method:
                self.unwrapped_phase = unwrap_qg(self.wrapped_phase,
                    self.module)
            else:
                print("Unrecognized option on unwrapping_method: %s" %
                    method)

            self.update_unwrapping_vis()
            self.update_overview_vis()
        else:
            print("Not unwrapping")


    @on_trait_change("unwrapping_vismode")
    def update_unwrapping_vis(self):
        if self.unwrapping_vismode == "phase":
            array = self.unwrapped_phase
            color = "spectral"
        elif self.unwrapping_vismode == "hibryd":
            array = self.unwrapped_phase
            color = "bone"
        else:
            print("Unrecognized option on unwrapping_vismode: %s" % 
                self.unwrapping_vismode)
            array = self.empty
            color = "prism"

        if self.plt_unwrapping is None:
            self.plt_unwrapping = self.scn_propagation.mlab.imshow(
                array, colormap=color,
                figure=self.scn_unwrapping.mayavi_scene)
        else:
            self.plt_unwrapping.mlab_source.set(scalars=array)


    ## PUT ALL-TOGHETER ##
    view = View(

        Group(
            HSplit(
                Group(
                    grp_datainput,
                    grp_ref_beam_parameters,
                    grp_mask_parameters,
                    grp_propagation_parameters,
                    grp_unwrapping_parameters,
                ),
                grp_overview_visualizer,
            ),
            label="Overview",
        ),

        Group(
            HSplit(
                Group(
                    grp_datainput,
                    grp_ref_beam_parameters,
                ),
                grp_ref_beam_visualizer,
            ),
            label="Reference beam",
        ),

        Group(
            HSplit(
                Group(
                    grp_datainput,
                    grp_mask_parameters,
                ),
                grp_mask_visualizer,
            ),
            label="Mask",
        ),

        Group(
            HSplit(
                Group(
                    grp_datainput,
                    grp_propagation_parameters,
                ),
                grp_propagation_visualizer,
            ),
            label="Propagation",
        ),


        Group(
            HSplit(
                Group(
                    grp_datainput,
                    grp_unwrapping_parameters,
                ),
                grp_unwrapping_visualizer,
            ),
            label="Unwrapping",
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
