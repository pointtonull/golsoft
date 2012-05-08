from traits.api import HasStrictTraits, Str, Int, Regex, List, Instance, Range
from traitsui.api import View, VGroup, Item, ListEditor

#-- Person Class ---------------------------------------------------------------

class Person(HasStrictTraits):

    # Trait definitions:
    name  = Str
    age   = Int
    phone = Regex(value = '000-0000', regex = '\d\d\d[-]\d\d\d\d')

    # Traits view definition:
    traits_view = View('name', 'age', 'phone',
                        width   = 0.18,
                        buttons = [ 'OK', 'Cancel' ])

#-- Sample Data ----------------------------------------------------------------

people = [
   Person(name = 'Dave Chomsky',        age = 39, phone = '555-1212'),
   Person(name = 'Mike Wakowski',       age = 28, phone = '555-3526'),
   Person(name = 'Joe Higginbotham',    age = 34, phone = '555-6943'),
   Person(name = 'Tom Derringer',       age = 22, phone = '555-7586'),
   Person(name = 'Dick Van Der Hooten', age = 63, phone = '555-3895'),
   Person(name = 'Harry McCallum',      age = 46, phone = '555-3285'),
   Person(name = 'Sally Johnson',       age = 43, phone = '555-8797'),
   Person(name = 'Fields Timberlawn',   age = 31, phone = '555-3547')
]

#-- ListEditorNotebookSelectionDemo Class --------------------------------------

class ListEditorNotebookSelectionDemo(HasStrictTraits):

    #-- Trait Definitions ------------------------------------------------------

    # List of people:
    people = List(Person)

    # The currently selected person:
    selected = Instance(Person)

    # The index of the currently selected person:
    index = Range(0, 7, mode = 'spinner')

    #-- Traits View Definitions ------------------------------------------------

    traits_view = View(
        Item('index'),
        '_',
        VGroup(
            Item('people@',
                  id         = 'notebook',
                  show_label = False,
                  editor     = ListEditor(use_notebook = True,
                                           deletable    = False,
                                           selected     = 'selected',
                                           export       = 'DockWindowShell',
                                           page_name    = '.name')
            )
        ),
        id   = 'traitsui.demo.Traits UI Demo.Advanced.'
               'List_editor_notebook_selection_demo',
        dock = 'horizontal')

    #-- Trait Event Handlers ---------------------------------------------------

    def _selected_changed(self, selected):
        self.index = self.people.index(selected)

    def _index_changed(self, index):
        self.selected = self.people[ index ]

#-- Set Up The Demo ------------------------------------------------------------

demo = ListEditorNotebookSelectionDemo(people = people)

if __name__ == "__main__":
    demo.configure_traits()
