#!/usr/bin/env python
#-*- coding: UTF-8 -*-


class Datum(property):
    def __init__(self, value):
        self.childs = set()
        self.value = value
        self.updated = True
        property.__init__(self, self.getter, self.setter)


    def set_property(self, updater, *dependencies):
        self.updater = updater

        for parent in dependencies:
            parent.childs.add(self)


    def deprecate(self):
        """
        Deprecates current value
        Propagates the deprecation to all childs
        """
        self.updated = False

        for child in self.childs:
            child.deprecate()


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
        Deprecates childs
        """

        if self.value != value:
            self.value = value

            for child in self.childs:
                child.deprecate()
        self.updated = True
        return value


class Depends:
    def __init__(self, *dependencies):
        """
        Decorator class
        """
        self.dependencies = dependencies
    
    def __call__(self, updater):
        datum = Datum(None)
        datum.deprecate()
        datum.set_property(updater, *self.dependencies)
        return datum



def main():

    class PEA(object):

        filename = Datum(0)
        @Depends(filename)
        def hologram(self):
            print("reading")
            return self.filename + 1
        
        @Depends(hologram)
        def spectrum(self):
            print("dft")
            return self.hologram + 1

        order_scale = Datum(0.8)
        zero_scale = Datum(1.2)
        @Depends(spectrum, order_scale, zero_scale)
        def masked_spectrum(self):
            print("enmasking")
            return self.spectrum + 1

        @Depends(spectrum)
        def distance(self):
            print("focusing")
            return self.spectrum + 1

        @Depends(distance)
        def propagation(self):
            print("creating propagation")
            return self.distance + 1

        @Depends(masked_spectrum, propagation)
        def propagated(self):
            print("propagating")
            return max(self.masked_spectrum, self.propagation) + 1

        @Depends(propagated)
        def reconstructed(self):
            print("idft")
            return self.propagated + 1

        @Depends(reconstructed)
        def module(self):
            print("module")
            return self.reconstructed + 1

        @Depends(reconstructed)
        def phase(self):
            print("phase")
            return self.reconstructed + 1

        unwrapper = Datum("quality guided")
        @Depends(phase, unwrapper)
        def unwrapped_phase(self):
            print("unwrapping")
            return self.phase + 1


    pea = PEA()
    print pea.hologram
    print pea.propagated
    pea.filename = 4
    print pea.propagated
    print pea.propagated
    print pea.hologram
    pea.filename = 3
    print pea.unwrapped_phase
    pea.order_scale = 1
    print pea.unwrapped_phase


if __name__ == "__main__":
    exit(main())
