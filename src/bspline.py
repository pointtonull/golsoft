#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from collections import namedtuple
from operator import itemgetter
import cairo
import numpy as np

tau = np.pi


class Point(namedtuple("Point", "x y")):
    __slots__ = ()

    def get_node(self, prev, post, softness=1./3):
        to_prev = self.distance(prev)
        to_post = self.distance(post)
        bisector = (prev * to_post + post * to_prev) / (to_prev + to_post)
        rel_bisector = (bisector - self)
        norm_bisector = rel_bisector / rel_bisector.distance()
        perp_normbisector = Point(norm_bisector.y, - norm_bisector.x)
        prev_handle = self + perp_normbisector * to_prev * softness
        post_handle = self - perp_normbisector * to_post * softness
        return prev_handle, self, post_handle

    def distance(self, another=None):
        another = another or Point(0, 0)
        diff = self - another
        dist = (diff.x ** 2 + diff.y ** 2) ** 0.5
        return dist

    def __sub__(self, another):
        return Point(self.x - another.x, self.y - another.y)

    def __add__(self, another):
        return Point(self.x + another.x, self.y + another.y)

    def __mul__(self, scalar):
        return Point(self.x * scalar, self.y * scalar)

    def __div__(self, scalar):
        return Point(self.x / scalar, self.y / scalar)


def polartorec(points, center=(256, 256)):
    xshift, yshift = center
    theta_scale = np.pi / xshift
    points = sorted(points, key=itemgetter(1), reverse=True)
    points = [
        (rho * np.cos(theta * theta_scale) + xshift,
         rho * np.sin(theta * theta_scale) + yshift)
        for rho, theta in points]
    return points


class Curve:

    def __init__(self, points, color=(1, 1, 1, 1)):
        self.points = [Point(x, y) for x, y in points]
        self.color = color

    def get_nodes(self):
        points = self.points
        nodes = [points[i+1].get_node(points[i], points[i+2])
            for i in range(-2, len(points) - 2)]
        return nodes



class Canvas:

    def __init__(self, shape=(512, 512), mode=cairo.FORMAT_RGB24):
        self.shape = shape
        self.mode = mode
        self.surface = cairo.ImageSurface(mode, shape[0], shape[1])
        self.context = cairo.Context(self.surface)
#        self.context.set_operator(cairo.OPERATOR_ADD)
        self.clear()


    def clear(self):
        oldoperator = self.context.get_operator()
        self.context.set_operator(cairo.OPERATOR_SOURCE)
        self.context.set_source_rgba(0, 0, 0, 0)
        self.context.paint()
        self.context.set_operator(oldoperator)


    def draw(self, curve, verbose=0):
        context = self.context
        nodes = curve.get_nodes()

        context.move_to(*curve.points[-2])
        for i in xrange(-1, len(nodes) - 1):
            control1 = nodes[i][-1]
            control2 = nodes[i + 1][0]
            final = nodes[i + 1][1]
            coords = [
                control1.x, control1.y,
                control2.x, control2.y,
                final.x, final.y,
            ]
            context.curve_to(*coords)
        context.set_source_rgba(*curve.color)
        context.fill()

        if verbose > 0:
            context.set_source_rgba(0, 0, 1, 1)
            for node in nodes:
                context.move_to(*node[0])
                context.line_to(*node[2])
                context.stroke()
            if verbose > 1:
                context.move_to(*curve.points[0])
                for point in curve.points:
                    context.line_to(*point)
                context.set_source_rgba(1, 0, 0, 1)
                context.close_path()
                context.stroke()

        return



    def as_array(self):
        data = self.surface.get_data()
        array = np.frombuffer(data, np.uint8)
        rows, cols = self.shape[:2]
        channels = array.size / (rows * cols)
        array.shape = (rows, cols, channels or None)
        return array


def main():
    import autopipe
    canvas = Canvas((512, 512))

    color = (.3, .3, .3, 1)
    curve = Curve(((75, 75, 100), (75, 325, 100), (325, 325, 100),
        (300, 75, 1/3.)), color)
    canvas.draw(curve, verbose=2)

    curve = Curve(((90, 90, 20), (90, 300, 100), (300, 200, 150)), color)
    canvas.draw(curve)

    curve = Curve(((120, 120, 110), (230, 260, 130)), color)
    canvas.draw(curve)

    array = canvas.as_array()
    autopipe.showimage(array)


if __name__ == "__main__":
    exit(main())
