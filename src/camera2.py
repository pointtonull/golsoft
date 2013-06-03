#!/usr/bin/env python
#-*- coding: UTF-8 -*-

import cv2.cv as cv
import cv2
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
    cap = cv.CaptureFromCAM(0)
    cv.NamedWindow("camera", cv.CV_WINDOW_NORMAL)
    cv.SetCaptureProperty(cap, cv.CV_CAP_PROP_FRAME_WIDTH, 720)
    cv.SetCaptureProperty(cap, cv.CV_CAP_PROP_FRAME_HEIGHT, 540)
    cols = int(cv.GetCaptureProperty(cap, cv.CV_CAP_PROP_FRAME_WIDTH))
    rows = int(cv.GetCaptureProperty(cap, cv.CV_CAP_PROP_FRAME_HEIGHT))
    grey = cv.CreateImage((cols, rows), 8, 1)
    cumulated = cv.CreateImage((cols, rows), 8, 1)

    equalize = True
    laplace = False

    settings = {
        "canny_avg": 10,
    }

    threshold1 = 600
    threshold2 = 200

    settings_names = sorted(settings.keys())
    setting_current = 0
    setting_name = settings_names[setting_current]

    while True:
        im = cv.QueryFrame(cap)
        cv.CvtColor(im, grey, cv.CV_BGR2GRAY)

        if equalize:
            cv.Smooth(grey, grey, param1=5, param2=5)
            cv.EqualizeHist(grey, grey)

        if laplace:
            cv.Canny(grey, grey, threshold1, threshold2)
            avg = cv.Avg(cumulated)[0]
            if avg > settings["canny_avg"] * 1.2:
                threshold1 *= 1.1
                threshold2 = threshold1 / 2.5
            if avg < settings["canny_avg"] / 1.2:
                threshold1 /= 1.1
                threshold2 = threshold1 / 2.5

        cv.ShowImage("camera", grey)

        key = cv.WaitKey(1)
        if key not in (-1, 1114085, 1245157): # None, block
            print("Key %d" % key)
            if key in ( # Capture one frame
                1048675, # c
                99, # c
                ):
                filenames = save_image(cap, 1)
                print("Capturing: %s" % ", ".join(list(filenames)))
            if key in ( # Capture ten frames
                1114179, # C
                1179715, # C (block)
                65603, # C
                131139, # C (block)
                ):
                filenames = save_image(cap, 10)
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
                settings[setting_name] += 1
                print("%s := %d" % (setting_name, settings[setting_name]))

            elif key in ( # Decrement value
                1113940, # Down
                65364,
                ):
                settings[setting_name] -= 1
                print("%s := %d" % (setting_name, settings[setting_name]))

            elif key in ( # Next setting
                1113939, # Right
                65363,
                ):
                setting_current = (setting_current + 1) % len(settings_names)
                setting_name = settings_names[setting_current]
                print("%s : %d" % (setting_name, settings[setting_name]))

            elif key in ( # Prev setting
                1113937, # Left
                65361,
                ):
                setting_current = (setting_current - 1) % len(settings_names)
                setting_name = settings_names[setting_current]
                print("%s : %d" % (setting_name, settings[setting_name]))

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
                break

if __name__ == "__main__":
    exit(main())
