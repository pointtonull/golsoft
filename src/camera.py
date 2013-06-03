#!/usr/bin/env python
#-*- coding: UTF-8 -*-

import time

from pygame import camera
import numpy as np
import pygame

from dft import get_shifted_idft, get_shifted_dft
from image import normalize, equalize
from automask import get_auto_mask
from propagation import  get_propagation_array



def get_spectrum(rhologram, **kw):
    spectrum = get_shifted_dft(rhologram)
    return spectrum


def apply_mask(spectrum, **kw):
    mask, masked, centered = get_auto_mask(spectrum, softness=0,
        radious_scale=1.5, zero_scale=1.6)
    return masked


def propagate(masked, **kw):
    propagation_array = get_propagation_array(masked.shape, kw["distance"])
    propagated = propagation_array * masked
    return propagated


def reconstruct(propagated, **kw):
    reconstructed = get_shifted_idft(propagated)
    return reconstructed


def get_complex_view(reconstructed, **kw):
    if kw["comp_view"] == "phase":
        return normalize(np.arctan2(reconstructed.real, reconstructed.imag))
    else:
        return np.abs(reconstructed)


def get_equalized(image, **kw):
    return equalize(image)


def get_normalized(image, **kw):
    return normalize(image)


def save_raw(cam):
    raw = cam.get_raw()
    filename = "cap-%d.raw" % time.time()
    file = open(filename, "w")
    file.write(raw)
    file.close()
    return filename


def main():
    pygame.init()
    camera.init()
    pygame.surfarray.use_arraytype("numpy")

    cams = camera.list_cameras()
    cam = camera.Camera(cams[0], (360, 296))
    cam = camera.Camera(cams[0], (640, 480))
    cam.start()
    fps = 25.0
    window = pygame.display.set_mode((640, 480), 0, 8)
    pygame.display.set_caption("Video")
    screen = pygame.display.get_surface()
    screen.set_palette([(i, i, i) for i in range(256)])

    print("Starting main loop")

    pea_list = [
        ("Spectrum", get_spectrum, get_equalized),
        ("Automask", apply_mask, get_normalized),
        ("Propagation", propagate, get_normalized),
        ("Reconstruction", reconstruct, get_complex_view),
        ]

    set_array = False
    set_equalize = False
    set_normalize = True
    set_pea = False
    pea_level = 1
    distance = 5
    comp_view = "phase"

    while True:
        events = pygame.event.get()
        for event in events:
            if event.type == pygame.QUIT:
                return
            elif event.type == pygame.KEYDOWN:
                if (event.key == pygame.K_q):
                    return

                # IMAGE PROCESSING
                elif (event.key == pygame.K_a):
                    set_array = not set_array
                    print("Converting to array: %s" % set_array)
                elif (event.key == pygame.K_n):
                    set_normalize = not set_normalize
                    print("Normalize: %s" % set_normalize)
                elif (event.key == pygame.K_e):
                    set_equalize = not set_equalize
                    print("Equalize: %s" % set_equalize)

                # PEA
                elif (event.key == pygame.K_p):
                    set_pea = not set_pea
                    print("PEA processing set: %s" % set_pea)
                    print("Setted pea to level %d, %s." % (
                        pea_level, pea_list[pea_level - 1][0]))
                elif (event.key == pygame.K_PAGEUP):
                    pea_level -= 1
                    pea_level = max(pea_level, 1)
                    print("Setted pea to level %d, %s." % (
                        pea_level, pea_list[pea_level - 1][0]))
                elif (event.key == pygame.K_PAGEDOWN):
                    pea_level += 1
                    pea_level = min(pea_level, len(pea_list))
                    print("Setted pea to level %d, %s." % (
                        pea_level, pea_list[pea_level - 1][0]))
                elif (event.key == pygame.K_TAB):
                    comp_view = "phase" if comp_view != "phase" else "mod"
                    print("PEA complex viewer set to: %s" % comp_view)

                # FOCUS DISTANCE
                elif (event.key == pygame.K_DOWN):
                    distance += 5
                    print("Distance: %.1f" % distance)
                elif (event.key == pygame.K_UP):
                    distance -= 5
                    print("Distance: %.1f" % distance)
                elif (event.key == pygame.K_LEFT):
                    distance -= .5
                    print("Distance: %.1f" % distance)
                elif (event.key == pygame.K_RIGHT):
                    distance += .5
                    print("Distance: %.1f" % distance)

                # FULSCREEN
                elif (event.key == pygame.K_f):
                    pygame.display.toggle_fullscreen()

                # CAPTURE
                elif (event.key == pygame.K_c):
                    filename = save_raw(cam)
                    print("Raw image saved to: %s" % filename)

        image = cam.get_image()

        if set_array:
            array = pygame.surfarray.array2d(image)

            if array.ndim > 2:
                array = round(array.mean(-1))
#                array = array[:,:,0] # red
#                array = array[:,:,0] # green
#                array = array[:,:,0] # blue

            if set_equalize:
                array = equalize(array).astype(int)
            elif set_normalize:
                array = normalize(array)

            pygame.surfarray.blit_array(screen, array)

        elif set_pea:
            array = pygame.surfarray.array2d(image)

            if array.ndim > 2:
                array = round(array.mean(-1))
#                array = array[:,:,0] # red
#                array = array[:,:,0] # green
#                array = array[:,:,0] # blue

            pea_algs = pea_list[:pea_level]
            pea_rep = pea_algs[-1][-1]

            for alg in pea_algs:
                try:
                    array = alg[1](array, distance=distance)
                except:
                    print("W: skipped framme's: %s" % alg[0])

            array = pea_rep(array, comp_view=comp_view).astype(int)

            pygame.surfarray.blit_array(screen, array)

        else:
            screen.blit(image, (0,0))

        pygame.display.flip()
        pygame.time.delay(int(1000./fps))

if __name__ == "__main__":
    exit(main())
