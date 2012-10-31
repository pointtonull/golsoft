#!/usr/bin/env python
#-*- coding: UTF-8 -*-

import cv2
import numpy as np
import time


def save_image(cap, cant=1):
    images = []
    # TODO: debe asegurarse de tomar %cant imágenes distintas
    for number in xrange(cant):
        image = cv.QueryFrame(cap)
        filename = "cap-%f.ppm" % time.time()
        images.append((filename, image))
        cv.SaveImage(filename, image)
        time.sleep(0.0)
        yield filename



def main():

    cam = cv2.VideoCapture(0)
    cv2.namedWindow("camera", 0) # Resizable
    readed, image = cam.read()

    grey = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    averaged = np.float32(grey)

    equalize = False
    laplace = False

    settings = {
        "canny_avg": 10,
        "acumulator_alpha": 0.1   
    }

    settings_names = sorted(settings.keys())
    setting_current = 0
    setting_name = settings_names[setting_current]

    threshold1 = 600.
    threshold2 = 200.

    while readed:

        readed, image = cam.read()
        cv2.cvtColor(image, cv2.COLOR_BGR2GRAY, grey)

        if equalize:
            cv2.equalizeHist(grey, grey)

        cv2.accumulateWeighted(grey, averaged, settings["acumulator_alpha"])
        normalized = cv2.convertScaleAbs(averaged)

        if laplace:
            normalized = cv2.Canny(normalized, threshold1, threshold2,
                normalized)
            avg = normalized.mean()
            scale = avg / settings["canny_avg"]
            if scale < 0.9 or scale > 1.1:
                scale = (scale + 4) / 5
                threshold1 *= scale
                threshold2 *= scale
                print("Canny := %f, %f" % (threshold1, threshold2))

        cv2.imshow('camera', normalized)

        key = cv2.waitKey(20)
        if key not in (-1, 1114085, 1245157): # None, block
            print("Key %d" % key)
            if key in ( # Capture one frame
                1048675, # c
                99, # c
                ):
                filenames = save_image(cam, 1)
                print("Capturing: %s" % ", ".join(list(filenames)))
            if key in ( # Capture ten frames
                1114179, # C
                1179715, # C (block)
                65603, # C
                131139, # C (block)
                ):
                filenames = save_image(cam, 10)
                print("Capturing: %s" % ", ".join(list(filenames)))

            elif key in ( # Toggle equalization
                1114181, # e
                1048677, # E
                1179717, # E (block)
                1245285, # e (block)
                101,     # e
                65605,   # E
                131141,  # E (block)
                196709,  # e (block)
                ):
                equalize = not equalize
                print("Equalize: %s" % equalize)

            elif key in ( # Toggle laplace
                1179724, # l
                1048684, # L (block(
                1114188, # L
                108, 
                65612,
                131148,
                196716,
                ):
                laplace = not laplace 
                print("Laplace: %s" % laplace)

            elif key in ( # Increment value
                1113938, # Up
                65362,
                ):
                settings[setting_name] *= 1.1
                print("%s := %f" % (setting_name, settings[setting_name]))

            elif key in ( # Decrement value
                1113940, # Down
                65364,
                ):
                settings[setting_name] /= 1.1
                print("%s := %f" % (setting_name, settings[setting_name]))

            elif key in ( # Next setting
                1113939, # Right
                65363,
                ):
                setting_current = (setting_current + 1) % len(settings_names)
                setting_name = settings_names[setting_current]
                print("%s : %f" % (setting_name, settings[setting_name]))

            elif key in ( # Prev setting
                1113937, # Left
                65361,
                ):
                setting_current = (setting_current - 1) % len(settings_names)
                setting_name = settings_names[setting_current]
                print("%s : %f" % (setting_name, settings[setting_name]))

            elif key in ( # Exit
                27, # ESC
                1048603, # ESC
                1114193, # q
                1048689, # Q
                1179729, # Q (block)
                1245297, # q (block)
                113,
                65617,
                131153,
                196721,
                ):
                readed = False

    cv2.destroyAllWindows()
    cam.release()

if __name__ == "__main__":
    exit(main())
