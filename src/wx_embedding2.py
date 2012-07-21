"""
This example shows to embed a Mayavi view in a wx frame.

The trick is to create a `HasTraits` object, as in the
mlab_traits_ui.py, mayavi_traits_ui.py, or the modifying_mlab_source.py
examples (:ref:`example_mlab_traits_ui`, :ref:`example_mayavi_traits_ui`,
:ref:`example_mlab_interactive_dialog`).

Calling the `edit_traits` method returns a `ui` object whose
`control` attribute is the wx widget. It can thus be embedded in a
standard wx application.

In this example, the wx part is very simple. See
:ref:`example_wx_mayavi_embed_in_notebook` for an example of more complex
embedding of Mayavi scenes in Wx applications.
"""

class MayaviView(HasTraits):

#-----------------------------------------------------------------------------
# Wx Code
import wx

class MainWindow(wx.Frame):

    def __init__(self, parent, id):
        wx.Frame.__init__(self, parent, id, 'Mayavi in Wx')
        self.mayavi_view = MayaviView()
        # Use traits to create a panel, and use it as the content of this
        # wx frame.
        self.control = self.mayavi_view.edit_traits(
                        parent=self,
                        kind='subpanel').control
        self.Show(True)

app = wx.PySimpleApp()
frame = MainWindow(None, wx.ID_ANY)
app.MainLoop()

