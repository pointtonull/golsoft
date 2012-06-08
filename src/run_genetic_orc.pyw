#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from itertools import izip

import numpy as np
from pyevolve import Crossovers, G2DList, GSimpleGA
from pyevolve import Mutators, Selectors, Initializators

from autopipe import showimage
from bspline import Canvas, Curve

from image import get_centered

pi = np.pi
tau = 2 * pi


def polartorec(points, center=(256, 256), scale=256):
    xshift, yshift = center

    rho_scale = 256 * 2 * points[-1][1]
    rho_scale = 256
    rhos = [rho * rho_scale for rho, theta, soft in points]

    thetas = [theta for rho, theta, soft in points]
    theta_scale = tau / sum(thetas)
    thetas = [sum(thetas[:i]) * theta_scale for i in range(len(thetas) - 1)]

    soft_scale = points[-1][-1] * 2
    soft_scale = 1
    softs = [soft for rho, theta, soft in points]

    points = [
        (
         rho * np.cos(theta) + xshift,
         rho * np.sin(theta) + yshift,
         soft * soft_scale
        )
        for rho, theta, soft in izip(rhos, thetas, softs)]
    return points[::-1]



class GeneticORC:
    def __init__(self, reference_image):
        self.reference = get_centered(reference_image)
        self.canvas = Canvas(reference_image.shape)
        self.rangemax = max(reference_image.shape)

        self.genome = G2DList.G2DList(10, 3) # rho, theta, softness
        self.genome.setParams(rangemin=0., rangemax=1., bestRawScore=0.00)

        self.genome.evaluator.set(self.fitness_func)
#        self.genome.evaluator.add(eval_func2)

        self.genome.initializator.set(Initializators.G2DListInitializatorReal)

#        self.genome.crossover.set(Crossovers.G2DListCrossoverSingleHPoint)
#        self.genome.crossover.set(Crossovers.G2DListCrossoverSingleVPoint)
        self.genome.crossover.set(Crossovers.G2DListCrossoverUniform)

#        self.genome.mutator.set(Mutators.G2DListMutatorIntegerRange)
        self.genome.mutator.set(Mutators.G2DListMutatorRealGaussian)
#        self.genome.mutator.set(Mutators.G2DListMutatorIntegerGaussian)

        self.best = None
        self.bestscore = None
        self.lastrun = None


    def run(self):
        darwin = GSimpleGA.GSimpleGA(self.genome)
#        darwin.selector.set(Selectors.GRankSelector)
#        darwin.selector.set(Selectors.GRouletteWheel)
        darwin.selector.set(Selectors.GTournamentSelector)
        darwin.setGenerations(100)
        darwin.setMutationRate(0.10)
#        darwin.terminationCriteria.set(GSimpleGA.ConvergenceCriteria)

        darwin.evolve(freq_stats=10)
        self.best = darwin.bestIndividual()
        self.lastrun = darwin
        return self.best


    def fitness_func(self, chromosome):
        # chromosome = [(rho, theta, softness), ...]
        self.canvas.clear()
        recpoints = polartorec(chromosome)
        curve = Curve(recpoints)
        self.canvas.draw(curve)
        array = self.canvas.as_array()
        array = array.mean(-1)
        array = get_centered(array)
        penalty = abs(array - self.reference) > 126
        penalty_sum = (penalty.sum() / self.rangemax) ** 2
        score = 10000. / penalty_sum
        if score > self.bestscore:
            self.bestscore = score
            self.canvas.clear()
            self.canvas.draw(curve, 1)
            array = self.canvas.as_array()
            print("New top score (%.2f / %.2f):" % (score, penalty_sum))
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
        print polartorec(best.genomeList)

    return 0

if __name__ == "__main__":
    exit(main())
