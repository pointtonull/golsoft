#!/usr/bin/env python
#-*- coding: UTF-8 -*-


"""
Toolbox to write hieralchical pipes. Is intented to realize the smallest
quantity of operations:

    * when a datum is updated all their dependents are deprecated

Datum is the base class that handle dependences, and propagate deprecation
    signals.
Depends is a helper decorator that enable easy creation of methods with
    Daturm ready interface.


The process must be a new-style class, e.g.:

    from dependences import Datum, Depends

    Class Salary(object):

        base = Datum(700)
        @Dependences(base)
        def retentions(self):
            return self.base * .2

        antiquity = Datum(5)
        @Dependences(base, antiquity)
        def gross(self):
            return self.base + self.base * .05 * self.antiquity

        childs = Datum(3)
        @Dependences(base, childs)
        def statements(self):
            return self.base * .15 + self.childs * 60

        @Dependences(gross, retentions, statements)
        def net(self):
            return self.gross - self.retentions + self.statements

"""


class Datum(property):
    def __init__(self, value=None):
        """
        Base class that allow the construction of complex hieralchical 
        proccess throgh the iteraction of simple elements and using a clear
        sintaxis.
        """
        self.dependences = set()
        self.dependents = set()
        self.value = value
        self.updated = value is not None
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
        self.value = value

        for dependent in self.dependents:
            dependent.deprecate()
        self.updated = True
        return value


class Depends:
    def __init__(self, *dependences):
        """
        Decorator class
        """
        self.dependences = dependences
    
    def __call__(self, updater):
        """
        Decorator trick that setups the datum updater and dependencies.
        """
        datum = Datum(None)
        datum.deprecate()
        datum.set_updater(updater)
        datum.set_dependences(*self.dependences)
        return datum

