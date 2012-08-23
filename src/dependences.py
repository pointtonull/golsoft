#!/usr/bin/env python
#-*- coding: UTF-8 -*-


class Datum(property):
    def __init__(self, value):
        """
        Base class that allow the construction of complex hieralchical 
        proccess throgh the iteraction of simple elements and using a clear
        sintaxis.
        """
        self.dependences = set()
        self.dependents = set()
        self.value = value
        self.updated = True
        property.__init__(self, self.getter, self.setter)


    def set_updater(self, updater):
        """
        Defines the function to be called as necessary to update the value.
        """
        
        assert callable(updater)
        self.updater = updater


    def set_dependences(self, *dependences):
        """
        Set the dependencies of this datum.
        Set this datum as dependent of its dependences.
        """

        for dependence in dependences:
            self.dependences.add(dependence)
            dependence.dependents.add(self)


    def deprecate(self):
        """
        Deprecates current value
        Propagates the deprecation to all dependents
        """
        self.updated = False

        for dependent in self.dependents:
            dependent.deprecate()


    def getter(self, parent):
        """
        If updated: returns value
        Else: updates and returns
        Else: raises ValueError
        """
        if self.updated:
            return self.value

        elif self.updater:
            return self.setter(parent, self.updater(parent))

        else:
            raise ValueError("No value nor updater.")            


    def setter(self, parent, value):
        """
        Uptades self value
        Deprecates dependents
        """
        if self.value != value:
            self.value = value

            for dependent in self.dependents:
                child.deprecate()
        self.updated = True
        return value


class Depends:
    def __init__(self, *dependences):
        """
        Decorator class
        """
        self.dependences = dependences
    
    def __call__(self, updater):
        datum = Datum(None)
        datum.deprecate()
        datum.set_property(updater, *self.dependences)
        return datum

