#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from autopipe import showimage
from bspline import Canvas, Curve, polartorec
from pyevolve import Crossovers, G2DList, GSimpleGA
from pyevolve import Mutators, Selectors
import numpy as np


from image import get_centered


class GeneticORC:
    def __init__(self, reference_image):
        self.reference = get_centered(reference_image)
        self.canvas = Canvas(reference_image.shape)
        self.rangemax = max(reference_image.shape)

        self.genome = G2DList.G2DList(10, 2)
        self.genome.setParams(rangemin=0, rangemax=self.rangemax,
            bestRawScore=0.00, roundDecimal=1)

        self.genome.evaluator.set(self.fitness_func)
#        self.genome.crossover.set(Crossovers.G2DListCrossoverSingleHPoint)
#        self.genome.crossover.set(Crossovers.G2DListCrossoverSingleVPoint)
        self.genome.crossover.set(Crossovers.G2DListCrossoverUniform)
#        self.genome.mutator.set(Mutators.G2DListMutatorIntegerRange)
        self.genome.mutator.set(Mutators.G2DListMutatorIntegerGaussian)

        self.best = None
        self.bestscore = None
        self.lastrun = None


    def run(self):
        darwin = GSimpleGA.GSimpleGA(self.genome)
        darwin.selector.set(Selectors.GRankSelector)
#        darwin.selector.set(Selectors.GRouletteWheel)
#        darwin.selector.set(Selectors.GTournamentSelector)
        darwin.setGenerations(100)
        darwin.setMutationRate(0.10)
#        darwin.terminationCriteria.set(GSimpleGA.ConvergenceCriteria)

        darwin.evolve(freq_stats=10)
        self.best = darwin.bestIndividual()
        self.lastrun = darwin
        return self.best


    def fitness_func(self, chromosome):
        #TODO para que el apareamiento mejore los ángulos deberían ser relativos
        self.canvas.clear()
        recpoints = polartorec(chromosome)
        curve = Curve(recpoints)
        points = curve.points
        distances = [points[i].distance(points[i+1])
            for i in xrange(-1, len(curve.points) - 1)]
        std_distances = (np.std(distances) / sum(distances)) ** 2
        self.canvas.draw(curve)
        array = self.canvas.as_array()
        array = array.mean(-1)
        array = get_centered(array)
        penalty = abs(array - self.reference) > 127
        penalty_sum = (penalty.sum() / self.rangemax) ** 2
        diff = abs(penalty_sum - std_distances)
        score = 1000. / (penalty_sum + std_distances + diff)
        if score > self.bestscore:
            self.bestscore = score
            self.canvas.clear()
            self.canvas.draw(curve, 2)
            array = self.canvas.as_array()
            print("New top score (%.2f / %.2f / %.2f):" % (score, std_distances,
                penalty_sum))
            showimage(array, penalty * 255)
        return score



def main():

    import sys
    from image import imread

    images = [(filename, imread(filename))
        for filename in sys.argv[1:]]

    for filename, image in images:
        print("File: %s" % filename)
        showimage(image)
        gorc = GeneticORC(image)
        best = gorc.run()
        print best
        print best.genomeList

    return 0

if __name__ == "__main__":
    exit(main())
