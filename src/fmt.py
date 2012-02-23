#!/usr/bin/env python
#-*- coding: UTF-8 -*-

"""
This is a simple implementation of the Fourier-Mellin Trasform
image = [m, n]
ft_magnitude = |fft(image)|
lp_ft_magnitude = logpolar(ft_magnitude)
fmt = fft(lp_ft_magnitude)
"""

from automask import get_mask
from autopipe import showimage
from enhance import equalize, logscale
from numpy import fft
from numpy import ndarray as atype
from numpy import sin, cos, exp, log
from scipy import misc
from scipy.ndimage import geometric_transform
import numpy as np

tau = 2 * np.pi


def get_shiftedfft(array):
    shiftedfft = fft.fftshift(fft.fft2(array))
    return shiftedfft


def get_logpolar(array, order=0):
    rows, cols = array.shape
    row0 = rows / 2.
    col0 = cols / 2.
    theta_scalar = tau / cols
    max_radius = (row0 ** 2 + col0 ** 2) ** .5
    rho_scalar = log(max_radius) / cols

    def out2in(dst_coords):
        theta, rho = dst_coords
        rho = exp(rho * rho_scalar)
        theta = theta * theta_scalar - tau / 2
        row_from = rho * cos(theta) + row0
        col_from = rho * sin(theta) + col0
        return row_from, col_from

    logpolar = geometric_transform(array, out2in, array.shape, order=order)
    return logpolar


def hi_pass_filter(array, radius=0.2):
    radius = min(array.shape) * radius
    window = np.bartlett(radius)
    mask = np.ones_like(array) - get_mask(array.shape, window)
    masked = array * mask
    return masked
    
#TODO: implementar cache de arrays YA
def get_fmt(array):
    """
    Follows this algoritm:
        * FFT with centered frecuencies
        * convolucionar la magnitud con un filtro high pass #TODO
        * Logpolar
        * FFT with centered frecuencies
    """
    fourier = get_shiftedfft(array)

    real_hi_passed = hi_pass_filter(fourier.real, .25)
    imag_hi_passed = hi_pass_filter(fourier.imag, .25)
    real_logpolar = get_logpolar(real_hi_passed, 3)
    imag_logpolar = get_logpolar(imag_hi_passed, 3)
    logpolar = real_logpolar + 1j * imag_logpolar
    fmt = get_shiftedfft(logpolar)
    return fmt


def get_correlation(image1, image2):
    """
    Todo esto es un invento y debe ser revisado
    """
    min_rows = min(image1.shape[0], image2.shape[0])
    min_cols = min(image1.shape[1], image2.shape[1])
    image1 = misc.imresize(image1, (min_rows, min_cols))
    image2 = misc.imresize(image2, (min_rows, min_cols))
    fmt1 = get_fmt(image1)
    fmt2 = get_fmt(image2)
    intensity1 = equalize(fmt1.real ** 2 + fmt1.imag ** 2)
    intensity2 = equalize(fmt2.real ** 2 + fmt2.imag ** 2)
    showimage(intensity1, intensity2)
    intensitydiff = (intensity2 - intensity1) ** 2
    showimage(equalize(intensitydiff))
    diff = intensitydiff.mean()
    correlation = (54**2 / (1 + diff)) ** 2
    return correlation


def get_fmt_correlation(image1, image2):
    fmt1 = get_fmt(image1)
    angle1 = np.angle(fmt1)
    fmt2 = get_fmt(image2)
    angle2 = np.angle(fmt2)
    theta_cross = exp((angle1 - angle2).imag)
    theta_phase = fft.ifft2(theta_cross).real
    theta_x, theta_y = np.unravel_index(theta_phase.argmax(), theta_phase.shape)
    dpp = 360 / theta_phase.shape[1]
    theta = dpp * (theta_y - 1)
    r1 = misc.imrotate(image1, -theta)
    r2 = misc.imrotate(image2, -(theta + 180))
    showimage(r1, r2)
    r1_f2 = get_shiftedfft(r1)
    angle1 = np.angle(get_shiftedfft(image1))
    angle2 = np.angle(r1_f2)
    r1_f2_cross = exp((angle1 - angle2).imag)
    r1_f2_phase = fft.ifft2(r1_f2_cross).real
    r2_f2 = get_shiftedfft(r2)    
    fft_image2 = get_shiftedfft(image2)
    angle1 = np.angle(fft_image2)
    angle2 = np.angle(r2_f2)
    r2_f2_cross = exp((angle1 - angle2).imag)
    r2_f2_phase = fft.ifft2(r2_f2_cross).real

    max_r1_f2 = r1_f2_phase.max()
    max_r2_f2 = r2_f2_phase.max()
    
    if max_r1_f2 > max_r2_f2:
        x, y = np.unravel_index(r1_f2_phase.argmax(), r1_f2_phase.shape)
        r = r1
    else:
        x, y = np.unravel_index(r2_f2_phase.argmax(), r2_f2_phase.shape)
        if theta < 180:
            theta = theta + 180
        else:
            theta = theta - 180
        r = r2
    print(r, x, y)
    
    tx = x - 1
    ty = y - 1
    
    if x > image1.shape[0] / 2:
        tx = tx - image1.shape[0]
    
    if y > image1.shape[1] / 2:
        ty = ty - image1.shape[1]
    
    input2_rectified = r
    move_ht = ty
    move_wd = tx

#    total_height = max(size(I1,1),(abs(move_ht)+size(input2_rectified,1)));
#    total_width =  max(size(I1,2),(abs(move_wd)+size(input2_rectified,2)));
#    combImage = zeros(total_height,total_width); registered1 = zeros(total_height,total_width); registered2 = zeros(total_height,total_width);

#    % if move_ht and move_wd are both POSITIVE
#    if((move_ht>=0)&&(move_wd>=0))
#        registered1(1:size(I1,1),1:size(I1,2)) = I1;
#        registered2((1+move_ht):(move_ht+size(input2_rectified,1)),(1+move_wd):(move_wd+size(input2_rectified,2))) = input2_rectified; 
#    elseif ((move_ht<0)&&(move_wd<0))   % if translations are both NEGATIVE
#        registered2(1:size(input2_rectified,1),1:size(input2_rectified,2)) = input2_rectified;
#        registered1((1+abs(move_ht)):(abs(move_ht)+size(I1,1)),(1+abs(move_wd)):(abs(move_wd)+size(I1,2))) = I1;
#    elseif ((move_ht>=0)&&(move_wd<0))
#        registered2((move_ht+1):(move_ht+size(input2_rectified,1)),1:size(input2_rectified,2)) = input2_rectified;
#        registered1(1:size(I1,1),(abs(move_wd)+1):(abs(move_wd)+size(I1,2))) = I1;
#    elseif ((move_ht<0)&&(move_wd>=0))
#        registered1((abs(move_ht)+1):(abs(move_ht)+size(I1,1)),1:size(I1,2)) = I1;
#        registered2(1:size(input2_rectified,1),(move_wd+1):(move_wd+size(input2_rectified,2))) = input2_rectified;    
#    end

#    if sum(sum(registered1==0)) > sum(sum(registered2==0))   % find the image with the greater number of zeros - we shall plant that one and then bleed in the other for the combined image
#        plant = registered1;    bleed = registered2;
#    else
#        plant = registered2;    bleed = registered1;
#    end

#    combImage = plant;
#    for p=1:total_height
#        for q=1:total_width
#            if (combImage(p,q)==0)
#                combImage(p,q) = bleed(p,q);
#            end
#        end
#    end
        
    
    
    
#    % Show final image
    
#    imshow(combImage, [0 255]);
