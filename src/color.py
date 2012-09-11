#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from autofocus import generic_minimizer
from image import normalize


def wavelength2rgb(wavelength, gamma=1., intensitymax=255):
    """
    Convert a wavelength on a rgb color.
    """

    wavelength = float(wavelength)

    if wavelength >= 350 and wavelength <= 439:
        red= -(wavelength - 440) / (440 - 350)
        green = 0.0
        blue = 1.0
    elif wavelength >= 440 and wavelength <= 489:
        red= 0
        green = (wavelength - 440) / (490 - 440)
        blue= 1.0
    elif wavelength >= 490 and wavelength <= 509:
        red = 0.0
        green = 1.0
        blue = -(wavelength - 510) / (510 - 490)
    elif wavelength >= 510 and wavelength <= 579:
        red = (wavelength - 510) / (580 - 510)
        green = 1.0
        blue = 0.0
    elif wavelength >= 580 and wavelength <= 644:
        red = 1.0
        green = -(wavelength - 645) / (645 - 580)
        blue = 0.0
    elif wavelength >= 645 and wavelength <= 780:
        red = 1.0
        green = 0.0
        blue = 0.0
    else:
        red = 0.0
        green = 0.0
        blue = 0.0

    if wavelength >= 350 and wavelength <= 419:
        factor = 0.3 + 0.7 * (wavelength - 350) / (420 - 350)
    elif wavelength >= 420 and wavelength <= 700:
         factor = 1.0
    elif wavelength >= 701 and wavelength <= 780:
        factor = 0.3 + 0.7 * (780 - wavelength) / (780 - 700)
    else:
        factor = 0.0
    
    red = adjustfactor(red, factor)
    green = adjustfactor(green , factor)
    blue = adjustfactor(blue, factor)

    return red, green, blue



def adjustfactor(color, factor, intensitymax=255, gamma=1.):
    if color == 0.:
        return 0
    else:
        return int(round(intensitymax * (color * factor) ** gamma))



def sature_color(rgbtuple):
    """
    A simple, not so acurate, way to normalize a dominant color.
    """
    red, green, blue = rgbtuple
    min_value = min(rgbtuple)
    max_value = max(rgbtuple)
    if min_value == max_value:
        return 0, 0, 0
    else:
        value_range = 255. / (max_value - min_value)
        red = int(round((red - min_value) * value_range))
        green = int(round((green - min_value) * value_range))
        blue = int(round((blue - min_value) * value_range))
        return red, green, blue



def rgb2wavelength(rgbtuple):
    """
    A regression wavelength2rgb reverter function.
    """
    original_tone = sature_color(rgbtuple)
    ored, ogreen, oblue = original_tone

    def distance(wavelength):
        newtone = sature_color(wavelength2rgb(wavelength))
        nred, ngreen, nblue = newtone
        distance = ((ored - nred) ** 2
            + (ogreen - ngreen) ** 2
            + (oblue - nblue) ** 2) ** .5
        return distance

    initial_guess = 580
    best_guess = generic_minimizer(distance, initial_guess)
    return best_guess



def guess_wavelength(image):
    """
    Helper function that try to deduce the image dominant wavelength.
    """

    if image.ndim == 2:
        print("W: Guess wavelength: a monocrome bitmap has not a dominant"
            " wavelength.")
        rgbcolor = (128, 128, 128)
    else:
        red = image[:, :, 0].mean()
        green = image[:, :, 1].mean()
        blue = image[:, :, 2].mean()
        rgbcolor = (red, green, blue)

    wavelength = rgb2wavelength(rgbcolor)
    return wavelength



def main():
    for wavelength in range(350, 800, 10):
        rgb = wavelength2rgb(wavelength)
        saturated = sature_color(rgb)
        nrgb = rgb2wavelength(rgb)
        print(wavelength, rgb, saturated, nrgb)

if __name__ == "__main__":
    exit(main())
