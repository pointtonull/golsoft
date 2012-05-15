from traits.api import HasStrictTraits, Str, Int, Regex, List, Instance
from traits.api import on_trait_change
from traits.api import Range
from traitsui.api import View, VGroup, Item, ListEditor


class Filter(HasStrictTraits):
    name  = Str
    age = Range(0, 100, mode='spinner')
    phone = Regex(value='000-0000', regex='\d\d\d[-]\d\d\d\d')

    traits_view = View(
        'name',
        'age',
        'phone',
        width = 0.18,
        buttons = ['OK', 'Cancel']
    )

pipe = [
    Filter(
        name = 'FFT',
        age = 41,
        phone = '555-1212'
    ),
    Filter(
        name = 'xRefBeam',
        age = 63,
        phone = '555-3895'
    ),
    Filter(
        name = 'Automask',
        age = 63,
        phone = '555-3895'
    ),
    Filter(
        name = 'Autofocus',
        age = 31,
        phone = '555-3547'
    ),
    Filter(
        name = 'Propagate',
        age = 46,
        phone = '555-3285'
    ),
    Filter(
        name = 'IDFT',
        age = 34,
        phone = '555-6943'
    ),
]


class ListEditorNotebookSelectionDemo(HasStrictTraits):
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
