#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from itertools import izip

import numpy as np
from pyevolve import Crossovers, G2DList, GSimpleGA
from pyevolve import Mutators, Selectors, Initializators, Consts

from autopipe import showimage
from bspline import Canvas, Curve
from image import get_centered, get_polar, normalize
from pea import generic_minimizer

pi = np.pi
tau = 2 * pi


class GeneticORC:
    def __init__(self, reference_image):
        self.reference = normalize(get_centered(reference_image))
        self.profile = (get_polar(self.reference) > 126).sum(-1)

        self.canvas = Canvas(reference_image.shape)
        self.rangemax = max(reference_image.shape)

        self.genome = G2DList.G2DList(10, 3) # rho, theta, softness
        self.genome.setParams(rangemin=0., rangemax=1., bestRawScore=0.00)

        self.genome.initializator.set(Initializators.G2DListInitializatorReal)

#        self.genome.mutator.set(Mutators.G2DListMutatorIntegerRange)
        self.genome.mutator.set(Mutators.G2DListMutatorRealGaussian)
#        self.genome.mutator.set(Mutators.G2DListMutatorIntegerGaussian)

        self.genome.evaluator.set(self.fitness_func)
#        self.genome.evaluator.add(eval_func2)

#        self.genome.crossover.set(Crossovers.G2DListCrossoverSingleHPoint)
#        self.genome.crossover.set(Crossovers.G2DListCrossoverSingleVPoint)
        self.genome.crossover.set(Crossovers.G2DListCrossoverUniform)
#        self.genome.minimax = Consts.minimaxType["minimize"]

        self.best = None
        self.bestscore = None
        self.lastrun = None


    def run(self):
        darwin = GSimpleGA.GSimpleGA(self.genome)
        darwin.selector.set(Selectors.GRankSelector)
#        darwin.selector.set(Selectors.GRouletteWheel)
#        darwin.selector.set(Selectors.GTournamentSelector)
        darwin.setGenerations(10)
        darwin.setMutationRate(0.10)
#        darwin.terminationCriteria.set(GSimpleGA.ConvergenceCriteria)

        darwin.evolve(freq_stats=10)
        self.best = darwin.bestIndividual()
        self.lastrun = darwin
        return self.best


    def buildpoints(self, points, center=(256, 256), scale=256):
        xshift, yshift = center

        thetas = [theta for rho, theta, soft in points]
        theta_scale = tau / sum(thetas)
        thetas = [sum(thetas[:i]) * theta_scale for i in range(len(thetas) - 1)]

        rhos = [rho for rho, theta, soft in points][:-1]
        rhos = [self.profile[int(round(thetas[i] * 512 / tau))]
            for i in xrange(len(thetas) - 1)]

        soft_scale = points[-1][-1] * 2
        softs = [soft for rho, theta, soft in points][:-1]

        points = [
            (
             rho * np.cos(theta) + xshift,
             rho * np.sin(theta) + yshift,
             soft * soft_scale
            )
            for rho, theta, soft in izip(rhos, thetas, softs)]

        return points[::-1]


    def fitness_func(self, chromosome):
        # chromosome = [(relrho, theta, softness), ...]
        points = self.buildpoints(chromosome)
        curve = Curve(points)

        self.canvas.clear()
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



class AdaptativeORC:
    def __init__(self, reference_image):
        self.reference = normalize(get_centered(reference_image))
        self.profile = (get_polar(self.reference) > 126).sum(-1)

        self.shape = reference_image.shape
        self.canvas = Canvas(self.shape)
        self.size = float(reference_image.size)

        self.best = None
        self.best_cost = None


    def buildpoints(self, points, center=(256, 256), scale=256):
        # chromosome = [(relrho, theta, softness), ...]
        xshift, yshift = center

        rows = points.size / 3
        points = points.reshape((rows, 3))
        thetas = [abs(theta) for rho, theta, soft in points]
        theta_scale = tau / sum(thetas)
        thetas = [sum(thetas[:i]) * theta_scale for i in range(len(thetas) - 1)]

        rhos = [rho * self.size ** 0.5 for rho, theta, soft in points][:-1]

#        try:
#            rhos = [self.profile[int(round(thetas[i] * 512 / tau))]
#                for i in xrange(len(thetas) - 1)]
#        except:
#            print(thetas, i)

        softs = [soft for rho, theta, soft in points][:-1]

        soft_scale = 1 + (points[-1][-1] * 2 - 1) ** 5

        points = [
            (
             rho * np.cos(theta) + xshift,
             rho * np.sin(theta) + yshift,
             soft * soft_scale
            )
            for rho, theta, soft in izip(rhos, thetas, softs)]

        return points[::-1]


    def cost_func(self, theta):
        # chromosome = [(relrho, theta, softness), ...]
        points = self.buildpoints(theta)
        curve = Curve(points)

        self.canvas.clear()
        self.canvas.draw(curve)
        array = self.canvas.as_array()
        array = array.mean(-1) # color -> grey tones
        array = get_centered(array)
        penalty = abs(array - self.reference) > 126
        cost = (penalty.sum() / self.size) * 1000
        cost += theta.std(0).sum()

        if cost < self.best_cost or self.best_cost is None:
            self.best_cost = cost
            self.best = theta
            self.canvas.clear()
            self.canvas.draw(curve, 1)
            array = self.canvas.as_array()
            print("New top score %.4f:" % cost)
            showimage(array, penalty * 255)

        return cost

    def run(self, initial_guest=None):
        initial_guest = initial_guest or np.random.random((7, 3))
        best = generic_minimizer(self.cost_func, initial_guest)
        return best
        


def main():

    import sys
    from image import imread

    images = [(filename, imread(filename))
        for filename in sys.argv[1:]]

    for filename, image in images:
        print("File: %s" % filename)
        showimage(image)
#        gorc = GeneticORC(image)
        gorc = AdaptativeORC(image)
        best = gorc.run()
        print best
        print gorc.buildpoints(best.genomeList)

    return 0

if __name__ == "__main__":
    exit(main())
