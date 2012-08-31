#!/usr/bin/env python
#-*- coding: UTF-8 -*-

import cv2.cv as cv
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

    equalize = False

    while True:
        im = cv.QueryFrame(cap)
        cv.CvtColor(im, grey, cv.CV_BGR2GRAY)

        if equalize:
            cv.EqualizeHist(grey, grey)


        cv.ShowImage("camera", grey)

        key = cv.WaitKey(1)
        if key not in (-1, 1114085, 1245157): # None, block
            print("Key %d" % key)
            if key in ( # Capture one frame
                1048675, # c
                ):
                filenames = save_image(cap, 1)
                print("Capturing: %s" % ", ".join(list(filenames)))
            if key in ( # Capture ten frames
                1114179, # C
                1179715, # C (block)
                ):
                filenames = save_image(cap, 10)
                print("Capturing: %s" % ", ".join(list(filenames)))
            elif key in ( # Toggle equalization
                1114181, # e
                1048677, # E
                1179717, # E (block)
                1245285, # e (block)
                ):
                equalize = not equalize
                print("Equialize: %s" % equalize)
            elif key in ( # Exit
                27, # ESC
                1048603, # ESC
                1114193, # q
                1048689, # Q
                1179729, # Q (block)
                1245297, # q (block)
                    ):
                break

if __name__ == "__main__":
    exit(main())
