#!/usr/bin/env python
#-*- coding: UTF-8 -*-

import numpy as np
from scipy import fftpack
import cv2

pi = np.pi
tau = 2 * pi

"""
                                 BRANCH CUTTING
                                 ==============
"""

def placecut(array, a, b, c, d, xsize, ysize, code):
    """
    Place a branch cut in the bitflags array from pixel (a,b) to pixel (c,d).
    The bit for the branch cut pixels is given by the value of "code".
    """

    # residue location is upper-left corner of 4-square
    if (c > a and a > 0):
        a += 1
    elif (c < a and c > 0):
        c += 1
    if (d > b and b > 0):
        b += 1
    elif (d < b and d > 0):
        d += 1

    if a == c and b == d:
      array[b * xsize + a] |= code
      return
    m = c - a if a < c else a - c
    n = d - b if b < d else b - d
    if (m > n):
      istep += 1 if a < c else -1
      r = float(d - b) / (c - a)
      for i in range(a, c + istep, istep):
        j = b + (i - a) * r + 0.5
        array[j * xsize + i] |= code
    else:
      jstep += 1 if b < d else -1
      r = float(c - a) / (d - b)
      for j in range(b, d + jstep, jstep):
        i = a + (j - b) * r + 0.5
        array[j * xsize + i] |= code
    return array


def disttoborder(bitflags, border_code, a, b, ra, rb, xsize, ysize):

    """
    Return the squared distance between the pixel (a,b) and the nearest border
    pixel.  The border pixels are encoded in the bitflags array by the value of
    "border_code".
    """

    ra = rb = 0
    for bs in range(xsize + ysize):
        found = False
        best_dist2 = 1000000 # just a large value
        # search boxes of increasing size until border pixel found
        for j in range(b - bs, b + bs):
            for i in range(a - bs, a + bs):
                k = j * xsize + i
                if (i ==0 or i >= xsize - 1 or j == 0 or j >= ysize - 1
                    or bitflags[k] & border_code):
                    found = True
                    dist2 = (j - b) * (j - b) + (i - a) * (i - a)
                    if dist2 < best_dist2:
                        best_dist2 = dist2
                        besta = i
                        bestb = j
        if found:
            ra = besta
            rb = bestb
            break
    return best_dist2



"""
                        FINDING AND ELIMINATING DIPOLES
                        ===============================
"""


def dipole(bitflags, xsize, ysize, branchcut_code):
    """
    Find the dipoles (adjoining residues of opposite sign) in the bitflags
    array, remove them, and place branch cut. The bits for the positive and
    negative residues are defined by POS_RES and NEG_RES (defined in brcut.h),
    while the branch cut bits are defined by "branchcut_code".
    """

    print("Dipoles...")
    for j in range(ysize):
        for i in range(xsize):
            k = j * xsize + i
            kk = 0
            if bitflags[k] & POS_RES:
                if i < xsize - 1 and bitflags[k+1] & NEG_RES:
                    kk = k + 1
                elif j < ysize-1:
                    if bitflags[k + xsize] & NEG_RES:
                        kk = k + xsize
            elif bitflags[k] & NEG_RES:
                if i < xsize - 1 and bitflags[k + 1] & POS_RES:
                    kk = k + 1
                elif j < ysize - 1:
                    if bitflags[k + xsize] & POS_RES:
                        kk = k + xsize
            if kk:
                print("Connecting dipoles %d,%d - %d,%d" % (i, j, kk % xsize,
                    kk / xsize))
                PlaceCut(bitflags, i, j, kk%xsize, kk/xsize, xsize, ysize,
                    branchcut_code)
                bitflags[k] &= ~RESIDUE
                bitflags[kk] &= ~RESIDUE


"""
                         EXTRACT PHASE FROM INPUT DATA
                         =============================
"""

def getphase(in_format, ifp, infile, phase, xsize, ysize):
    """
    Get the phase data from the input data.
    Allocate the temporary memory required and read the input file.
    in_format =
        0 for 8-byte complex data,
        1 for 4-byte complex data,
        2 for 1-byte phase data,
        3 for 4-byte float phase
    data Output: normalized (0..1) phase values in array "phase"
    """

    print("Reading input data...")
    if in_format == 1:   # short int real and imaginary parts
        ReadShort(ifp, in_data, 2 * xsize * ysize, infile)
    elif in_format == 0: # floating-pt real and imag parts
        ReadFloat(ifp, in_data, 2 * xsize * ysize, infile)
    elif in_format == 3: # floating-pt phase
        ReadFloat(ifp, in_data, xsize * ysize, infile)
    else:                # 1-byte quantized phase
        ReadByte(ifp, in_data, xsize * ysize, infile)

    return extractphase(in_format, in_data, phase, xsize, ysize, 1)


def extractphase(in_format, in_data, phase, xsize, ysize, status_flag):
    """
    Extract the phase from the input data.  If status_flag is 0 then do not
    print any status.
    """

    if in_format==0 or in_format == 1:
        for j in range(ysize):
            if status_flag and j % 100 == 0:
                print("%d " % j)
                if j + 100 >= ysize:
                    print("")
            for i in range(xsize):
                if in_format == 1:
                    x = in4_data[2 * (j * xsize + i)]
                    y = in4_data[2 * (j * xsize + i) + 1]
                else:
                    x = in8_data[2 * (j * xsize + i)]
                    y = in8_data[2 * (j * xsize + i) + 1]
                r = (x**2 + y**2) ** .5
                r = 0. if r == 0.0 else 1. / r


"""
                        MAKE SURFACE CONGRUENT TO PHASE 
                        ===============================
"""

def congruentsoln(surf, phase, qual_map, xsize, ysize):
    """
    Make the surface in the array "surf" congruent to the wrapped
    phase (scaled to 0-1) in the array "phase".  Ignore the      
    zero-weights in the quality map array qual_map.              
    """
    # ensure surface is positive
    rmin = surf[0]
    for k in xrange(xsize * ysize):
        if rmin > surf[k]:
            rmin = surf[k]

    npos = -rmin + 2.0 if rmin < 0 else 0.0
    for x in xrange(xsize * ysize):
        surf[k] += npos

    # find best additive offset (to minimize discontinuities)
    for kk in xrange(11):
        rr = 0.1 * kk
        # create temporary surface by adding a constant rr
        for k in xrange(xsize * ysize):
            if qual_map and qual_map[k] == 0.0:
                continue # ignore 0-wts
            n = surf[k] + rr
            r = surf[k] + rr - n
            r -= phase[k]
            if r < -0.5:
                ss[k] = n - 1 + phase[k]
            elif r > 0.5:
                ss[k] = n + 1 + phase[k]
            else:
                ss[k] = n + phase[k]

        # measure the discontinuities in the temporary surface
        rnum = 0
        for k in xrange(xsize * ysize):
            if qual_map and qual_map[k] == 0.0:
                continue # ignore 0-wts
            i = k % xsize
            j = k / xsize
            if i > 0:
                if qual_map and qual_map[k-1] == 0.0:
                    continue # ignore
                r = ss[k] - ss[k-1]
                if r > 0.5 or r < -0.5:
                    rnum += 1

            if j > 0:
                if qual_map and qual_map[k - xsize] == 0.0:
                    continue # ignore
                r = ss[k] - ss[k-xsize]
                if r > 0.5 or r < -0.5:
                    rnum += 1

        # save the rr that gives the minimum error (rnum)
        if kk == 0 or rnum < rnum_min:
            rnum_min = rnum
            rr_min = rr

        print("Offset %f yields %d discontinuities" % (rr, rnum))

    delete(ss)
    print("Adding offset %lf to surface" % rr_min)
    for k in xrange(xsize * ysize):
        surf[k] += rr_min

    for k in xrange(xsize * ysize):
        n = surf[k]
        r = surf[k] - n
        r -= phase[k]
        if r < -0.5:
            surf[k] = n - 1 + phase[k]
        elif r > 0.5:
            surf[k] = n + 1 + phase[k]
        else:
            surf[k] = n + phase[k]
        surf[k] -= npos # also un-do addition of npos


"""
              FUNCTIONS FOR COMPUTING PHASE DERIVATIVES IN ARRAYS
              ===================================================
"""

def dxphasegradient(phase, dx, xsize, ysize):
    """
    Compute the wrapped row differences (d/dx partial derivatives of the
    phase) and store them in the array "dx".
    """

    for j in xrange(ysize):
        for i in xrange(xsize):
            k = j * xsize + i
            if i < xsize - 1: 
                grad = gradient(phase[k + 1], phase[k])
            else:
                grad = Gradient(phase[k - 1], phase[k])
            dx[k] = grad


def dyphasegradient(phase, dy, xsize, ysize):
    """
    Compute the wrapped column differences (d/dy partial derivatives) and
    store them in the array "dy".
    """
    for j in xrange(ysize):
        for i in xrange(xsize):
            k = j * xsize + i
            if j < ysize - 1:
                grad = Gradient(phase[k + xsize], phase[k])
            else:
                grad = Gradient(phase[k-xsize], phase[k])
            dy[k] = grad


"""
                         EXTRACT PHASE FROM INPUT DATA
                         =============================
"""

def getphase(filename):
    """
    Get the phase data from the input data.    Allocate the temporary
    memory required and read the input file.

    in_format = 0 for 8-byte complex data,
        1 for 4-byte complex data,
        2 for 1-byte phase data,
        3 for 4-byte float phase data

    Output: normalized (0 - 1) phase values in array "phase" 
    """
    image = imread(filename)
    phase = extractphase(image)
    return phase


def extractphase(image, verbose=0):
    """
    Extract the phase from the input data. If status_flag is 0
    then do not print any status.
    """
    #TODO: no usar o reemplazar con implementación en PEA
    return


"""
               FUNCTIONS FOR fLYNN'S MIN. DISCONTINUITY ALGORITHM
               ==================================================
"""

"""

Main iteration (and initialization steps) for Flynn's minimum discontinuity
phase-unwrapping algorithm.    The quality values (weights) are stored in the
qual_map array.    This array and the phase array are xsize x ysize pixels in
size.    The bitflags array, which stores the change and direction info, and
the other arrays (value, vjump and hjump) are working arrays for managing the
data.  The dimensions of these arrays are (xsize + 1) x (ysize + 1).  If there
are masked (border) pixels, they should be coded as zero weights in qual_map.

void FlynnMinDisc(float *phase, int *value, unsigned char *bitflags,
                                    float *qual_map, short *vjump, short *hjump,
                                    int xsize, int ysize)
{
    double                    eps = 1.0e-06, phd
    short                     jmpold
    int                         j, i, jnext, inext
    int                         ilast, jlast
    int                         new_loops, new_edges
    int                         val, value_incr
    int                         loop_found
    int                         value_change

    # Compute phase jump counts.
     * If dx(i,j) and dy(i,j) represent the wrapped derivatives
     * phase(i+1,j) - phase(i,j) and phase(i,j+1) - phase(i,j),
     * respectively, then the jump count hjump(i,j) corresponds to
     * dy(i-1,j-1), and vjump(i,j) corresponds to dx(i-1,j-1).

    for (j=1; j <= ysize-1; j += 1) {
        for (i=1; i <= xsize; i += 1) {
            phd = phase[j*xsize + i-1] - phase[(j-1)*xsize + i-1] + eps
            hjump[j*(xsize+1) + i] = (short) nint(phd)


    for (j=1; j <= ysize; j += 1) {
        for (i=1; i <= xsize-1; i += 1) {
            phd = phase[(j-1)*xsize + i] - phase[(j-1)*xsize + i-1] + eps
            vjump[j*(xsize+1) + i] = (short) nint(phd)


    # Make add nodes bitflags initially
    for (j=0; j <= ysize; j += 1) {
        for (i=0; i <= xsize; i += 1) {
            bitflags[j*(xsize+1) + i] |= THIS_TIME


    # Main iteration
    do {
        new_loops = 0
        new_edges = 0
        # Add left-to-right edges to tree
        for (j=0; j <= ysize-1; j += 1) {
            for (i=0; i <= xsize; i += 1) {
                jnext = j+1
                inext = i
                if ((bitflags[j*(xsize+1) + i] & THIS_TIME)
                            or (bitflags[jnext*(xsize+1) + inext] & THIS_TIME)) {
                    if (jnext <= 1 or jnext >= ysize 
                             or inext <= 0 or inext >= xsize) {
                        value_incr = -isign(1, -vjump[jnext*(xsize+1) + inext])

                    else {
                        # quality map values must be between 0 and 1
                        val = (int)(1.0 + BIG*qual_map[jnext*xsize + inext])
                        value_incr    
                             = -isign(val, -vjump[jnext*(xsize+1) + inext])

                    value_change = value[j*(xsize+1) + i] + (int)value_incr
                                                                     - value[jnext*(xsize+1) + inext]
                    if (value_change > 0) {
                        new_edges += 1
                        # revise values in subtree of [jnext][inext]
                        # and check for loop
                        ilast = i;    jlast = j;    loop_found = 0
                        ChangeExten(inext, jnext, ilast, jlast, &loop_found,
                                                value_change, value, bitflags, xsize, ysize)
                        if (loop_found) {
                            RemoveLoop(inext, jnext, i, j, value, bitflags,
                                                 vjump, hjump, xsize, ysize)
                            ChangeOrphan(inext, jnext, -value_change, value,
                                                     bitflags, xsize, ysize)
                            new_loops += 1

                        else {
                            # add edge and remove edges to [jnext][inext]
                            bitflags[j*(xsize+1) + i] |= RIGHT
                            if (jnext<ysize) 
                                bitflags[(jnext+1)*(xsize+1) + inext] &= ~LEFT
                            if (inext<xsize) 
                                bitflags[jnext*(xsize+1) + inext+1] &= ~UP
                            if (inext > 0)
                                bitflags[jnext*(xsize+1) + inext-1] &= ~DOWN





        # Add top-to-bottom edges to tree
        for (j=0; j <= ysize; j += 1) {
            for (i=0; i <= xsize-1; i += 1) {
                jnext = j
                inext = i+1
                if ((bitflags[j*(xsize+1) + i] & THIS_TIME)
                         or (bitflags[jnext*(xsize+1) + inext] & THIS_TIME)){
                    if (jnext <= 0 or jnext >= ysize 
                             or inext <= 1 or inext >= xsize) {
                        value_incr = -isign(1, hjump[jnext*(xsize+1) + inext])

                    else {
                        # quality map values must be between 0 and 1
                        val = (int)(1.0 + BIG*qual_map[jnext*xsize + inext])
                        value_incr 
                                = -isign(val, hjump[jnext*(xsize+1) + inext])

                    value_change = value[j*(xsize+1) + i] + (int)value_incr
                                                                 - value[jnext*(xsize+1) + inext]
                    if (value_change > 0) {
                        new_edges += 1
                        ilast = i; jlast = j; loop_found = 0
                        ChangeExten(inext, jnext, ilast, jlast, &loop_found,
                                                value_change, value, bitflags, xsize, ysize)
                        if (loop_found){
                            RemoveLoop(inext, jnext, i, j, value, bitflags,
                                                 vjump, hjump, xsize, ysize)
                            ChangeOrphan(inext, jnext, -value_change, value,
                                                     bitflags, xsize, ysize)
                            new_loops += 1

                        else {
                            bitflags[j*(xsize+1) + i] |= DOWN
                            if (jnext<ysize) 
                                bitflags[(jnext+1)*(xsize+1) + inext] &= ~LEFT
                            if (jnext > 0) 
                                bitflags[(jnext-1)*(xsize+1) + inext] &= ~RIGHT
                            if (inext<xsize) 
                                bitflags[jnext*(xsize+1) + inext+1] &= ~UP





        # Add right-to-left edges to tree
        for (j=ysize; j >=1 ; j--) {
            for (i=xsize; i >=0; i--) {
                jnext = j-1
                inext = i
                if ((bitflags[j*(xsize+1) + i] & THIS_TIME)
                             or (bitflags[jnext*(xsize+1) + inext] & THIS_TIME)){
                    if (j <= 1 or j >= ysize or i <= 0 or i >= xsize) {
                        value_incr = -isign(1, vjump[j*(xsize+1) + i])

                    else {
                        # quality map values must be between 0 and 1
                        val = (int)(1.0 + BIG*qual_map[j*xsize + i])
                        value_incr = -isign(val, vjump[j*(xsize+1) + i])

                    value_change = value[j*(xsize+1) + i] + (int)value_incr
                                                                 - value[jnext*(xsize+1) + inext]
                    if (value_change > 0) {
                        new_edges += 1
                        ilast = i; jlast = j; loop_found = 0
                        ChangeExten(inext, jnext, ilast, jlast, &loop_found,
                                                value_change, value, bitflags, xsize, ysize)
                        if (loop_found){
                            RemoveLoop(inext, jnext, i, j, value, bitflags, 
                                                 vjump, hjump, xsize, ysize)
                            ChangeOrphan(inext, jnext, -value_change, value,
                                                     bitflags, xsize, ysize)
                            new_loops += 1

                        else {
                            bitflags[j*(xsize+1) + i] |= LEFT
                            if (jnext > 0)
                                bitflags[(jnext-1)*(xsize+1) + inext] &= ~RIGHT
                            if (inext < xsize)
                                bitflags[jnext*(xsize+1) + inext+1] &= ~UP
                            if (inext > 0) 
                                bitflags[jnext*(xsize+1) + inext-1] &= ~DOWN





        # Add bottom-to-top edges to tree
        for (j=ysize; j >=0; j--) {
            for (i=xsize; i >= 1; i--) {
                jnext = j
                inext = i-1
                if ((bitflags[j*(xsize+1) + i] & THIS_TIME)
                             or (bitflags[jnext*(xsize+1) + inext] & THIS_TIME)){
                    if (j <= 0 or j >= ysize or i <= 1 or i >= xsize) {
                        value_incr = -isign(1, -hjump[j*(xsize+1) + i])

                    else {
                        # quality map values must be between 0 and 1
                        val = (int)(1.0 + BIG*qual_map[j*xsize + i])
                        value_incr = -isign(val, -hjump[j*(xsize+1) + i])

                    value_change = value[j*(xsize+1) + i] + (int)value_incr
                                                                 - value[jnext*(xsize+1) + inext]
                    if (value_change > 0) {
                        new_edges += 1
                        ilast = i; jlast = j; loop_found = 0
                        ChangeExten(inext, jnext, ilast, jlast, &loop_found,
                                                value_change, value, bitflags, xsize, ysize)
                        if (loop_found){
                            RemoveLoop(inext, jnext, i, j, value, bitflags,
                                                 vjump, hjump, xsize, ysize)
                            ChangeOrphan(inext, jnext, -value_change, value,
                                                     bitflags, xsize, ysize)
                            new_loops += 1

                        else {
                            bitflags[j*(xsize+1) + i] |= UP
                            if (jnext < ysize) 
                                bitflags[(jnext+1)*(xsize+1) + inext] &= ~LEFT
                            if (jnext > 0) 
                                bitflags[(jnext-1)*(xsize+1) + inext] &= ~RIGHT
                            if (inext > 0) 
                                bitflags[jnext*(xsize+1) + inext-1] &= ~DOWN





        print("New edges: %d    New loops: %d",new_edges,new_loops)
        # Update testability matrix
        for (j=0; j <= ysize; j += 1) {
            for (i=0; i <= xsize; i += 1) {
                if (bitflags[j*(xsize+1) + i] & NEXT_TIME) {
                    bitflags[j*(xsize+1) + i] |= THIS_TIME

                else {
                    bitflags[j*(xsize+1) + i] &= ~NEXT_TIME
                    bitflags[j*(xsize+1) + i] &= ~THIS_TIME



 while (new_edges > 0)

    # Compute new unwrapped phase and exit
    print("Computing revised unwrapped phase...")
    for (i = 0; i < xsize - 1; i += 1) {
        phd = phase[0*xsize + (i + 1)] - phase[0*xsize + i] + eps
        jmpold = (short) nint(phd)
        phase[0*xsize + (i + 1)] 
                += (float)(vjump[1*(xsize + 1) + (i + 1)] - jmpold)

    for (j = 0; j < ysize - 1; j += 1) {
        for (i = 0; i < xsize; i += 1) {
            phd = phase[(j + 1)*xsize + i] 
                                        - phase[j*xsize + i] + eps
            jmpold = (short) nint(phd)
            phase[(j + 1)*xsize + i] 
                += (float)(hjump[(j + 1)*(xsize + 1) + (i + 1)] - jmpold)



#
 *    fmg.c -- function for multigrid solution of weighted 
 *                     least-squares phase-unwrapping problem

#include <stdio.h>
#include <math.h>
#include "fmg.h"
#include "grid.h"
#include "dxdygrad.h"

#    Call the functions for performing the multigrid phase-
#    unwrapping algorithm.    If dywts is a null pointer, then
#    dxwts are copied into an array for dywts.
void MultigridUnwrap(float *soln, float *dx, float *dy, float *dxwts,
    float *dywts, int xsize, int ysize, int num_cycles, int num_iter)
{
    int    n, dywts_was_null=0, coarsest_dim=3
    if (dywts==NULL) {
        dywts_was_null = 1
        AllocateDouble(&dywts, xsize*ysize, "dy wts")
        for (n=0; n<xsize*ysize; n += 1) dywts[n] = dxwts[n]

    for (n=0; n<num_cycles; n += 1) {
        print("\nFMG CYCLE %d", n+1)
        FullMultigridVcycle(soln, dx, dy, dxwts, dywts,
                                                xsize, ysize, num_iter, coarsest_dim)

    if (dywts_was_null) delete(dywts)

#
 * getqual.c -- generate or process quality map and set quality mode

#include <stdio.h>
#include <math.h>
#include "getqual.h"
#include "file.h"
#include "util.h"
#include "pi.h"
#include "qualvar.h"
#include "qualpseu.h"
#include "qualgrad.h"

# Generate a new quality map from the phase data, or process
# the quality map that was input.
void GetQualityMap(int mode, float *qual_map, float *phase, 
                                     unsigned char *bitflags, int border_code,
                                     int tsize, int xsize, int ysize)
{
    float    *temp
    double rmin, rmax, rscale
    int        i, j, k
    # process phase gradients
    if (mode==variance) {
        AllocateFloat(&temp, xsize*ysize, "temp data")
        PhaseVariance(phase, qual_map, bitflags, border_code, temp,
                                    tsize, xsize, ysize)
        delete(temp)
        # convert from cost to quality, and scale to interval (0,1)
        for (rmin = rmax = qual_map[0], k=0; k<xsize*ysize; k += 1) {
            if (rmin > qual_map[k]) rmin = qual_map[k]
            if (rmax < qual_map[k]) rmax = qual_map[k]

        print("Min & max phase derivative variance = %lf, %lf",
                     rmin, rmax)
        rscale = (rmin != rmax) ? 1.0/(rmax - rmin) : 0.0
        for (k=0; k<xsize*ysize; k += 1) {
            qual_map[k] = (rmax - qual_map[k])*rscale
            if (bitflags and (bitflags[k]&border_code)) qual_map[k] = 0.0


    else if (mode==gradient) {
        AllocateFloat(&temp, xsize*ysize, "temp data")
        MaxPhaseGradients(phase, qual_map, bitflags, border_code, temp,
                                            tsize, xsize, ysize)
        delete(temp)
        # convert from cost to quality, and scale to interval (0,1)
        for (rmin = rmax = qual_map[0], k=0; k<xsize*ysize; k += 1) {
            if (rmin > qual_map[k]) rmin = qual_map[k]
            if (rmax < qual_map[k]) rmax = qual_map[k]

        print("Min&max 'max phase gradient' = %lf %lf", rmin,rmax)
        rscale = (rmin != rmax) ? 1.0/(rmax - rmin) : 0.0
        for (k=0; k<xsize*ysize; k += 1) {
            qual_map[k] = (rmax - qual_map[k])*rscale
            if (bitflags and (bitflags[k]&border_code)) qual_map[k] = 0.0


    else if (mode==pseudocorr) {
        AllocateFloat(&temp, xsize*ysize, "temp data")
        PseudoCorrelation(phase, qual_map, bitflags, border_code, temp,
                                            tsize, xsize, ysize)
        for (k=0; k<xsize*ysize; k += 1) {
            if (bitflags and (bitflags[k]&border_code)) qual_map[k] = 0.0

        delete(temp)

    else if (mode==none) {
        for (k=0; k<xsize*ysize; k += 1) {
            qual_map[k] = 1.0
            if (bitflags and (bitflags[k]&border_code)) qual_map[k] = 0.0


    else { # quality map was input
        # scale to interval (0,1)
        for (rmin = rmax = qual_map[0], k=0; k<xsize*ysize; k += 1) {
            if (rmin > qual_map[k]) rmin = qual_map[k]
            if (rmax < qual_map[k]) rmax = qual_map[k]

        print("Min & max corr. coeff. = %lf, %lf", rmin, rmax)
        rscale = (rmin != rmax) ? 1.0/(rmax - rmin) : 0.0
        for (k=0; k<xsize*ysize; k += 1) {
            qual_map[k] = (qual_map[k] - rmin)*rscale
            if (bitflags and (bitflags[k]&border_code)) qual_map[k] = 0.0




# Determine the quality mode based on a keyword, and return a
# corresponding integer.    If keyword unrecognized, return -1.
# Allow the keyword "none" only if allow_none = 1.
int SetQualityMode(char *modekey, char *qualfile, int allow_none)
{
    int mode
    if (Keyword(modekey, "none")) {
        mode = none
        print("No weights supplied")
        if (!allow_none) {
            fprintf(stderr, "Invalid mode: %s", modekey)
            mode = -1;    # "none" not allowed


    else if (Keyword(modekey, "max_corr")) {
        mode = corr_coeffs
        if (Keyword(qualfile, "none")) {
            fprintf(stderr, "Must supply quality image for this mode")
            mode = -1


    else if (Keyword(modekey, "min_grad")) {
        if (!Keyword(qualfile, "none"))
            print("Correlation image file will be ignored")
        mode = gradient

    else if (Keyword(modekey, "min_var")) {
        if (!Keyword(qualfile, "none"))
            print("Correlation image file will be ignored")
        mode = variance

    else if (Keyword(modekey, "max_pseu")) {
        if (!Keyword(qualfile, "none"))
            print("Correlation image file will be ignored")
        mode = pseudocorr

    else {
        fprintf(stderr, "Invalid mode = '%s'.    Must be", modekey); 
        fprintf(stderr, "max_corr, min_var, max_grad, max_pseu "
                                        "or none.")
        mode = -1;    # invalid

    return mode
"""

"""
                 GENERATE BRANCH CUTS BY GOLDSTEIN'S ALGORITHM
                 =============================================
"""

def goldsteinbranchcuts(bitflags, maxcutlen, numres, xsize, ysize, branchcut_code):

    """
    Goldstein's phase-unwrapping algorithm. The bitflags store the masked pixels
    (to be ignored) and the residues and accumulates other info such as the branch
    cut pixels.
    """

    bench = ysize / 100.
    if bench < 1:
        bench = 1
    if MaxCutLen < 2:
        MaxCutLen = 2
    max_active = numres + 10
    active_list = max_active + 1 # book keeping data
    # branch cuts
    print("Computing branch cuts")
    for j in xrange(ysize):
        if j % bench == 0:
            print("%d ", j / bench)
        for i in xrange(xsize):
            k = j * xsize + i
            if (bitflags[k] & (POS_RES | NEG_RES)) and not(bitflags[k] & VISITED):
                bitflags[k] |= "VISITED" # turn on visited flag
                bitflags[k] |= "ACTIVE" # turn on active flag
                charge = 1 if (bitflags[k] & POS_RES) else -1
                num_active = 0
                num_active += 1
                active_list[num_active] = k
                if num_active > max_active:
                    num_active = max_active
                for boxsize in xrange(3, 2 * maxcutlen + 1, 2):
                    bs2 = boxsize / 2
                    for ka in xrange(num_active):
                        boxctr_i = active_list[ka] % xsize
                        boxctr_j = active_list[ka] / xsize
                        for jj in xrange(boxctr_j - bs2, boxctr_j + bs2 + 1):
                            for ii in xrange(boxctr_i - bs2, boxctr_i + bs2 + 1):
                                kk = jj * xsize + ii
                                if ii < 0 or ii >= xsize or jj < 0 or jj >= ysize:
                                    continue
                                else:
                                    if (ii == 0 or ii == xsize - 1 or jj == 0
                                        or jj == ysize - 1
                                        or (bitflags[kk] & BORDER)):
                                        charge = 0
                                        disttoborder(bitflags, BORDER, boxctr_i,
                                            boxctr_j, ri, rj, xsize, ysize)
                                        placecut(bitflags, ri, rj, boxctr_i,
                                            boxctr_j, xsize, ysize, branchcut_code)

                                    elif ((bitflags[kk] & (POS_RES | NEG_RES))
                                        and not(bitflags[kk] & ACTIVE)):
                                        if not(bitflags[kk] & VISITED):
                                            charge += 1 if bitflags[kk] & POS_RES else -1
                                            bitflags[kk] |= "VISITED" # set flag
                                        num_active += 1
                                        active_list[num_active] = kk
                                        if num_active > max_active:
                                            num_active = max_active
                                        bitflags[kk] |= ACTIVE # set active flag
                                        placecut(bitflags, ii, jj, boxctr_i,
                                            boxctr_j, xsize, ysize, branchcut_code)
                                    if charge == 0:
                                        # mark all active pixels inactive
                                        for ka in xrange(num_active): 
                                            bitflags[active_list[ka]] &= ~ACTIVE
                                            # turn flag off
                                        return
    # if (bitflags ...
    # else
     # for (ii ...
     # for (jj ...
    # for (ka ...
     # for (boxsize ...

                if (charge != 0): # connect branch cuts to rim
                    min_dist = xsize + ysize # large value
                    for ka in xrange(num_active):
                        ii = active_list[ka] % xsize
                        jj = active_list[ka] / xsize
                        if ((dist == gisttoborder(bitflags, BORDER,
                            ii, jj, ri, rj, xsize, ysize)) < min_dist):
                            min_dist = dist
                            near_i = ii
                            near_j = jj
                            rim_i = ri
                            rim_j = rj

                    placecut(bitflags, near_i, near_j, rim_i, rim_j, xsize,
                        ysize, branchcut_code)
 
                # mark all active pixels inactive
                for ka in xrange(num_active): 
                    bitflags[active_list[ka]] &= ~ACTIVE # turn flag off
    # if (bitflags ...
    # for (i ...
    # for (j ...
    print("")
    return


"""
       FUNCTION FOR COMPUTING PHASE DERIVATIVE (WRAPPED PHASE DIFFERENCE)
       ==================================================================


#include <stdio.h>
#include <math.h>
#include "grad.h"
# return wrapped phase difference
float Gradient(float p1, float p2)
{
    float    r
    r = p1 - p2
    if (r > 0.5) r -= 1.0
    if (r < -0.5) r += 1.0
    return r

#
 * grid.c - functions for (weighted) multigrid phase unwrapping 

#include <stdio.h>
#include <math.h>
#include "grid.h"
#include "gridops.h"
#include "gridmem.h"
#include "relax.h"

# Main function for weighted multigrid (called recursively)
void FullMultigridVcycle(float *soln, float *dx, float *dy,
     float *dxwts, float *dywts, int w, int h, int numit, int mindim)
{
    float    *dx2=NULL, *dy2=NULL, *soln2=NULL
    float    *dxwts2=NULL, *dywts2=NULL
    int        w2 = w/2, h2 = h/2
    if (!Coarsest(w, h, mindim)) {
        dxwts2 = Allocate(w2, h2, dxwts_type)
        dywts2 = Allocate(w2, h2, dywts_type)
        dx2 = Allocate(w2, h2, dx_type)
        dy2 = Allocate(w2, h2, dy_type)
        soln2 = Allocate(w2, h2, soln_type)
        RestrictDxwts(dxwts2, w2, h2, dxwts, w, h)
        RestrictDywts(dywts2, w2, h2, dywts, w, h)
        Restrict(dx2, dy2, w2, h2, dx, dy, dxwts, dywts, soln, w, h); 
        Zero(soln2, w2, h2)
        FullMultigridVcycle(soln2, dx2, dy2, dxwts2, dywts2, w2, h2,
                                                numit, mindim)
        ProlongAndAccumulate(soln, w, h, soln2, w2, h2, dxwts2, dywts2)

    # perform V-cycle multigrid on fine grid
    Vcycle(soln, dx, dy, dxwts, dywts, w, h, numit, mindim)


# V-cycle multigrid algorithm (called recursively)
void Vcycle(float *soln, float *dx, float *dy, float *dxwts,
                        float *dywts, int w, int h, int numit, int mindim)
{
    float *dx2=NULL, *dy2=NULL, *soln2=NULL
    float *dxwts2=NULL, *dywts2=NULL
    int        w2 = w/2, h2 = h/2
    if (!Coarsest(w, h, mindim)) {
        Relax(soln, dx, dy, dxwts, dywts, w, h, numit)
        dxwts2 = Allocate(w2, h2, dxwts_type)
        dywts2 = Allocate(w2, h2, dywts_type)
        dx2 = Allocate(w2, h2, dx_type)
        dy2 = Allocate(w2, h2, dy_type)
        soln2 = Allocate(w2, h2, soln_type)
        RestrictDxwts(dxwts2, w2, h2, dxwts, w, h)
        RestrictDywts(dywts2, w2, h2, dywts, w, h)
        Restrict(dx2, dy2, w2, h2, dx, dy, dxwts, dywts, soln, w, h); 
        Zero(soln2, w2, h2)
        Vcycle(soln2, dx2, dy2, dxwts2, dywts2, w2, h2, numit, mindim); 
        ProlongAndAccumulate(soln, w, h, soln2, w2, h2, dxwts2, dywts2)
        Relax(soln, dx, dy, dxwts, dywts, w, h, numit)

    else { # coarsest
        Relax(soln, dx, dy, dxwts, dywts, w, h, 2*w*h); 



# multigrid restriction operation for dx weights
void RestrictDxwts(float *coarse, int wc, int hc, 
                                     float *fine, int wf, int hf) 
{
    int         a, b, i, j, k, m, n
    int         k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k20
    float     *f, *c
    f = fine;                 c = coarse
    #    Indexes:
    #    k1     k2        k3
    #    k4     k5=k    k6
    #    k7     k8        k9
    for (j=0; j<hc; j += 1) {
        b = 2*j
        for (i=0; i<wc; i += 1) {
            a = 2*i
            k = b*wf + a
            k5 = k
            k6 = (a < wf - 1) ? k5 + 1 : k5 - 1
            k4 = (a > 0) ? k5 - 1 : k5 + 1
            if (b < hf - 1) {
                k7 = k4 + wf;     k8 = k5 + wf;        k9 = k6 + wf

            else {
                k7 = k4 - wf;     k8 = k5 - wf;        k9 = k6 - wf

            if (b > 0) {
                k1 = k4 - wf;     k2 = k5 - wf;     k3 = k6 - wf

            else {
                k1 = k4 + wf;     k2 = k5 + wf;     k3 = k6 + wf

            if (f[k5]==0.0 and f[k9]==0.0) c[j*wc + i] = 0.0
            else if (f[k6]==0.0 and f[k8]==0.0) c[j*wc + i] = 0.0
            else if (f[k6]==0.0 and f[k9]==0.0) c[j*wc + i] = 0.0
            else if (f[k5]==0.0 and f[k8]==0.0) c[j*wc + i] = 0.0
            else if ((i==0 or c[j*wc + i - 1] != 0.0) &&
                ((f[k4]==0.0 and f[k8]==0.0) 
                    or    (f[k5]==0.0 and f[k7]==0.0)))
                                                                         c[j*wc + i] = 0.0
            else c[j*wc + i] = 0.25*(f[k5] + f[k6] + f[k8] + f[k9])




# multigrid restriction operation for dy weights
void RestrictDywts(float *coarse, int wc, int hc, 
                                     float *fine, int wf, int hf) 
{
    int         a, b, i, j, k, m, n
    int         k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k20
    float    *f, *c
    f = fine;                 c = coarse
    #    Indexes:
    #    k1     k2        k3
    #    k4     k5=k    k6
    #    k7     k8        k9
    for (j=0; j<hc; j += 1) {
        b = 2*j
        for (i=0; i<wc; i += 1) {
            a = 2*i
            k = b*wf + a
            k5 = k
            k6 = (a < wf - 1) ? k5 + 1 : k5 - 1
            k4 = (a > 0) ? k5 - 1 : k5 + 1
            if (b < hf - 1) {
                k7 = k4 + wf;        k8 = k5 + wf;        k9 = k6 + wf

            else {
                k7 = k4 - wf;        k8 = k5 - wf;        k9 = k6 - wf

            if (b > 0) {
                k1 = k4 - wf;        k2 = k5 - wf;        k3 = k6 - wf

            else {
                k1 = k4 + wf;        k2 = k5 + wf;        k3 = k6 + wf

            if (f[k5]==0.0 and f[k9]==0.0) c[j*wc + i] = 0.0
            else if (f[k6]==0.0 and f[k8]==0.0) c[j*wc + i] = 0.0
            else if (f[k8]==0.0 and f[k9]==0.0) c[j*wc + i] = 0.0
            else if (f[k5]==0.0 and f[k6]==0.0) c[j*wc + i] = 0.0
            else if ((j==0 or c[(j - 1)*wc + i] != 0.0) &&
                    ((f[k2]==0.0 and f[k6]==0.0)
                        or    (f[k3]==0.0 and f[k5]==0.0)))
                                                                         c[j*wc + i] = 0.0
            else c[j*wc + i] = 0.25*(f[k5] + f[k6] + f[k8] + f[k9])


         
#
 * gridmem.c - dynamic memory management functions
 *                         for multigrid phase unwrapping 

#include <stdio.h>
#include <math.h>
#include "gridmem.h"

static int sizes[NUM_TYPES][MAX_ARRAYS]
static int num[NUM_TYPES]
static float *arrays[NUM_TYPES][MAX_ARRAYS]

# allocate array: if already allocated, return pointer to it
float *Allocate(int w, int h, ArrayType type)
{
    int size = w*h, i, j, k
    k = (int) type
    for (i=0; i<num[k]; i += 1) {
        if (sizes[k][i]==size) 
            return arrays[k][i]

    sizes[k][num[k]] = size
    arrays[k][num[k]] = (float *) malloc(size*sizeof(float))
    if (!arrays[k][num[k]]) {
        print("Error: cannot allocate memory (%d bytes)", size*4)
        return 0;    # error

    return arrays[k][num[k] += 1]


# free all allocated arrays
void FreeAll()
{ 
    int i, j
    for (j=0; j<NUM_TYPES; j += 1) {
        for (i=0; i<num[j]; i += 1) {
            delete(arrays[j][i])



#
 * gridops.c - prolongation and restriction operators for 
 *                         multigrid phase unwrapping.    Also contains
 *                         functions for detecting coarsest array size
 *                         and initializing an array to zero.

#include <stdio.h>
#include <math.h>
#include "gridops.h"
#define MIN(x,y)    (((x) < (y)) ? (x) : (y))

# Multigrid prolongation operator: weighted or unweighted.
# Resamples coarse grid values and ADDS to fine grid.
# (Note: for unweighted case, pass null pointers for
# coarse_dxwts and coarse_dywts.)
void ProlongAndAccumulate(float *fine, int wf, int hf,
                                                    float *coarse, int wc, int hc,
                                                    float *coarse_dxwts, float *coarse_dywts)
{
    int        i, j, k, k1, k2, k3, k4, a, b, c
    float    x1, x2, x3, x4, y1, y2, y3, y4, w1, w2, w3, w4, w
    float    w1x, w2x, w3x, w4x, w1y, w2y, w3y, w4y

    for (j=0, k=0; j<hc; j += 1) {
        for (i=0; i<wc; i += 1, k += 1) {
            a = 2*i
            b = 2*j
            c = b*wf + a
            k1 = k
            k2 = (i < wc - 1) ? k1 + 1 : k1
            k3 = (j < hc - 1) ? k1 + wc : k1
            k4 = (i < wc - 1) ? ((j < hc - 1) ? k1 + wc + 1 : k2) : k3
            x1 = coarse[k1]
            x2 = coarse[k2]
            x3 = coarse[k3]
            x4 = coarse[k4]
            if (coarse_dxwts and coarse_dywts) {
                w1x = coarse_dxwts[k1]
                w2x = coarse_dxwts[k2]
                w3x = coarse_dxwts[k3]
                w4x = coarse_dxwts[k4]
                w1y = coarse_dywts[k1]
                w2y = coarse_dywts[k2]
                w3y = coarse_dywts[k3]
                w4y = coarse_dywts[k4]
                w1 = MIN(w1x, w1y)
                w2 = MIN(w2x, w2y)
                w3 = MIN(w3x, w3y)
                w4 = MIN(w4x, w4y)
                y1 = x1
                w = w1 + w2
                y2 = (w > 0.0) ? (w1*x1 + w2*x2)/w : 0.5*(x1 + x2)
                w = w1 + w3
                y3 = (w > 0.0) ? (w1*x1 + w3*x3)/w : 0.5*(x1 + x3)
                w = w1 + w2 + w3 + w4
                y4 = (w > 0.0) ? (w1*x1 + w2*x2 + w3*x3 + w4*x4)/w
                                                                    : 0.25*(x1 + x2 + x3 + x4)

            else {
                y1 = x1;                                     y2 = 0.5*(x1 + x2)
                y3 = 0.5*(x1 + x3);                y4 = 0.25*(x1 + x2 + x3 + x4)

            fine[c] += y1
            fine[c + 1] += y2
            fine[c + wf] += y3
            fine[c + wf + 1] += y4
            # what if wf is odd? fill in extra columns
            if (i==wc - 1) {
                b = 2*j
                for (a=2*wc; a<wf; a += 1) {
                    c = b*wf + a
                    fine[c] += y2
                    fine[c + wf] += y4


            # what if hf is odd? fill in extra rows
            if (j==hc - 1) {
                a = 2*i
                for (b=2*hc; b<hf; b += 1) {
                    c = b*wf + a
                    fine[c] += y3
                    fine[c + 1] += y4


            # fill in lower right-hand corner
            if (i==wc - 1 and j==hc - 1) {
                for (b=2*hc; b<hf; b += 1) {
                    for (a=2*wc; a<wf; a += 1) {
                        c = b*wf + a
                        fine[c] += y4







# Multigrid restriction operator: weighted or unweighted.
# Restricts residuals of derivatives to coarser grid.
# (Note: for unweighted case, pass null pointers for
# dxwts and dywts.)
void Restrict(float *dx2, float *dy2, int wc, int hc,
                            float *dx, float *dy, float *dxwts, float *dywts,
                            float *soln, int wf, int hf) 
{
    int         a, b, i, j, k
    int         k1, k2, k3, k4, k5, k6, k7, k8, k9
    int         k1x, k2x, k3x, k4x, k5x, k6x, k7x, k8x, k9x
    int         k1y, k2y, k3y, k4y, k5y, k6y, k7y, k8y, k9y
    float     scale
    double    f1, f2, f3, f4, f5, f6, f7, f8, f9, wmult
    double    w1, w2, w3, w4, w5, w6, w7, w8, w9
    w1 = w2 = w3 = w4 = w5 = w6 = w7 = w8 = w9 = 1.0

    #    DX
    scale = (wf - 1.0)/(wc - 1.0)
    for (j=0; j<hc; j += 1) {
        b = 2*j
        for (i=0; i<wc-1; i += 1) {     # note: a < wf - 1
            a = 2*i
            k = b*wf + a
            #    Indexes:
            #    k1    k2    k3
            #    k4    k5    k6
            #    k7    k8    k9
            k5 = k
            k4 = (a > 0) ? k5 - 1 : k5 + 1
            k6 = (a < wf - 1) ? k5 + 1 : k5 - 1
            k2 = (b > 0) ? k5 - wf : k5 + wf
            k8 = (b < hf - 1) ? k5 + wf : k5 - wf
            k1 = (a > 0) ? k2 - 1 : k2 + 1
            k3 = (a < wf - 1) ? k2 + 1 : k2 - 1
            k7 = (a > 0) ? k8 - 1 : k8 + 1
            k9 = (a < wf - 1) ? k8 + 1 : k8 - 1
            k1x = k2;                 k4x = k5;                 k7x = k8
            k2x = k3;                 k5x = k6;                 k8x = k9
            k3x = k3 + 1;         k6x = k6 + 1;         k9x = k9 + 1
            if (dxwts and dywts) {
                # weights: compute as MIN(k, k-1)
                if ((w1 = dxwts[k1x]) > dxwts[k1]) w1 = dxwts[k1]
                if ((w2 = dxwts[k2x]) > dxwts[k2]) w2 = dxwts[k2]
                if ((w3 = dxwts[k3x]) > dxwts[k3]) w3 = dxwts[k3]
                if ((w4 = dxwts[k4x]) > dxwts[k4]) w4 = dxwts[k4]
                if ((w5 = dxwts[k5x]) > dxwts[k5]) w5 = dxwts[k5]
                if ((w6 = dxwts[k6x]) > dxwts[k6]) w6 = dxwts[k6]
                if ((w7 = dxwts[k7x]) > dxwts[k7]) w7 = dxwts[k7]
                if ((w8 = dxwts[k8x]) > dxwts[k8]) w8 = dxwts[k8]
                if ((w9 = dxwts[k9x]) > dxwts[k9]) w9 = dxwts[k9]
                wmult = 0.25*w5 + 0.125*(w2 + w8 + w4 + w6) 
                                                    + 0.0625*(w1 + w3 + w7 + w9)

            # dx residuals
            f1 = w1*(dx[k1] - (soln[k1x] - soln[k1])); 
            f2 = w2*(dx[k2] - (soln[k2x] - soln[k2])); 
            f3 = w3*(dx[k3] - (soln[k3x] - soln[k3])); 
            f4 = w4*(dx[k4] - (soln[k4x] - soln[k4])); 
            f5 = w5*(dx[k5] - (soln[k5x] - soln[k5])); 
            f6 = w6*(dx[k6] - (soln[k6x] - soln[k6])); 
            f7 = w7*(dx[k7] - (soln[k7x] - soln[k7])); 
            f8 = w8*(dx[k8] - (soln[k8x] - soln[k8])); 
            f9 = w9*(dx[k9] - (soln[k9x] - soln[k9])); 
            if (dxwts and dywts) {
                if (wmult > 1.0e-6) dx2[j*wc + i] = scale*(0.25*f5 
                                 +    0.125*(f4 + f6 + f2 + f8)
                                                +    0.0625*(f1 + f3 + f7 + f9))/wmult
                else dx2[j*wc + i] = dx[k]

            else {
                dx2[j*wc + i] = scale*(0.25*f5 + 0.125*(f4 + f6 + f2 + f8)
                                                                     +    0.0625*(f1 + f3 + f7 + f9))



    # correct dx at boundary
    for (j=0; j<hc; j += 1) {
        dx2[j*wc + (wc - 1)] = -dx2[j*wc + (wc - 2)]


    #    DY
    scale = (hf - 1.0)/(hc - 1.0)
    for (j=1; j<hc; j += 1) {     # note: b > 0
        b = 2*j
        for (i=0; i<wc; i += 1) {
            a = 2*i
            k = b*wf + a
            #    Indexes:
            #    k1    k2    k3
            #    k4    k5    k6
            #    k7    k8    k9
            k5 = k
            k4 = (a > 0) ? k5 - 1 : k5 + 1
            k6 = (a < wf - 1) ? k5 + 1 : k5 - 1
            k2 = (b > 0) ? k5 - wf : k5 + wf
            k8 = (b < hf - 1) ? k5 + wf : k5 - wf
            k1 = (a > 0) ? k2 - 1 : k2 + 1
            k3 = (a < wf - 1) ? k2 + 1 : k2 - 1
            k7 = (a > 0) ? k8 - 1 : k8 + 1
            k9 = (a < wf - 1) ? k8 + 1 : k8 - 1
            k1y = k4;                 k2y = k5;                 k3y = k6; 
            k4y = k7;                 k5y = k8;                 k6y = k9
            k7y = k7 + wf;        k8y = k8 + wf;        k9y = k9 + wf
            if (dxwts and dywts) {
                # weights: compute as MIN(k, k-1)
                if ((w1 = dywts[k1y]) > dywts[k1]) w1 = dywts[k1]
                if ((w2 = dywts[k2y]) > dywts[k2]) w2 = dywts[k2]
                if ((w3 = dywts[k3y]) > dywts[k3]) w3 = dywts[k3]
                if ((w4 = dywts[k4y]) > dywts[k4]) w4 = dywts[k4]
                if ((w5 = dywts[k5y]) > dywts[k5]) w5 = dywts[k5]
                if ((w6 = dywts[k6y]) > dywts[k6]) w6 = dywts[k6]
                if ((w7 = dywts[k7y]) > dywts[k7]) w7 = dywts[k7]
                if ((w8 = dywts[k8y]) > dywts[k8]) w8 = dywts[k8]
                if ((w9 = dywts[k9y]) > dywts[k9]) w9 = dywts[k9]
                wmult = 0.25*w5 + 0.125*(w2 + w8 + w4 + w6) 
                                                    + 0.0625*(w1 + w3 + w7 + w9)

            # dy residuals
            f1 = w1*(dy[k1] - (soln[k1y] - soln[k1])); 
            f2 = w2*(dy[k2] - (soln[k2y] - soln[k2])); 
            f3 = w3*(dy[k3] - (soln[k3y] - soln[k3])); 
            f4 = w4*(dy[k4] - (soln[k4y] - soln[k4])); 
            f5 = w5*(dy[k5] - (soln[k5y] - soln[k5])); 
            f6 = w6*(dy[k6] - (soln[k6y] - soln[k6])); 
            f7 = w7*(dy[k7] - (soln[k7y] - soln[k7])); 
            f8 = w8*(dy[k8] - (soln[k8y] - soln[k8])); 
            f9 = w9*(dy[k9] - (soln[k9y] - soln[k9])); 
            if (dxwts and dywts) {
                if (wmult > 1.0e-6) dy2[j*wc + i] = scale*(0.25*f5 
                                 +    0.125*(f4 + f6 + f2 + f8)
                                                +    0.0625*(f1 + f3 + f7 + f9))/wmult
                else dy2[j*wc + i] = scale*dy[k]

            else {
                dy2[j*wc + i] = scale*(0.25*f5 + 0.125*(f4 + f6 + f2 + f8)
                                                                     +    0.0625*(f1 + f3 + f7 + f9))



    # correct dy at boundary
    for (i=0; i<wc; i += 1) {
        dy2[(hc - 1)*wc + i] = -dy2[(hc - 2)*wc + i]



# returns 1 if w or h is less than 3
int Coarsest(int w, int h, int mindim) 
{
    if (mindim==0)
        return (w/2 < 3 or h/2 < 3)
    else
        return (w/2 < mindim or h/2 < mindim)


# initialize array to zeros
void Zero(float *soln, int w, int h)
{
    int n=w*h, k
    for (k=0; k<n; k += 1) {
        soln[k] = 0.0


#
 *    histo.c -- functions for histogramming and thresholding

#include <stdio.h>
#include <math.h>
#include "histo.h"
#define NUM_BINS                     10000
#define NUM_BINS_TO_PRINT    10

#
 * Compute histogram, apply histogram equalization and
 * apply threshold.    If thresh=0, the threshold is determined
 * automatically.    If percent_flag is 1, the threshold is
 * selected to threshold the percentage of the pixels
 * defined by "percentage".    The bitflags array is "ignore",
 * and ignore_code defines those pixels to ignore in the
 * array "qual" of values to be histogrammed and thresholded.
 * The results overwrite the values in "qual".

void HistoAndThresh(float *qual, int xsize, int ysize,
                             double thresh, int percent_flag, double percentage, 
                             unsigned char *ignore, int ignore_code)
{
    int         i, j, k
    int         num_zero=0, sum, num
    int         histo[NUM_BINS + 1]
    double    first_val, last_val, rmult
    print("Computing histogram, equalizing & applying threshold")
    # compute and print histogram
    Histogram(qual, xsize, ysize, histo, NUM_BINS,
                        ignore, ignore_code)
    # stretch histogram so at least 5% in first & last bins
    for (num=0.0, i=0; i<NUM_BINS; i += 1) { # find number of pixels
        num += histo[i]

    for (sum=0.0, i=0; i<NUM_BINS; i += 1) {
        sum += histo[i]
        if (sum >= 0.05*num) break

    first_val = (1.0*i)/NUM_BINS
    for (sum=0.0, i=0; i<NUM_BINS; i += 1) {
        sum += histo[i]
        if (sum >= 0.95*num) break

    last_val = (1.0*i)/NUM_BINS
    if (last_val <= first_val) last_val = first_val + 1.0/NUM_BINS
    rmult = 1.0/(last_val - first_val)
    # stretch weights
    for (j=0; j<xsize*ysize; j += 1) {
        if (ignore[j] & ignore_code) continue;     # ignore these
        if (qual[j] < first_val)
            qual[j] = 0.0
        else if (qual[j] > last_val)
            qual[j] = 1.0
        else
            qual[j] = (qual[j] - first_val)*rmult

    # compute histogram again
    print("Stretching histogram...")
    Histogram(qual, xsize, ysize, histo,
                        NUM_BINS, ignore, ignore_code)
    for (num=0.0, i=0; i<NUM_BINS; i += 1) {    # find number of pixels
        num += histo[i]

    if (thresh==0.0) {
        # Find optimal threshold by finding the trough of the
        # histogram.    Failing that, choose percentage = 50.
        int    downhill=0, uphill=0, num_coarse_bins = 10
        double    coarse_histo[NUM_BINS]
        print("Looking for trough of histogram to set threshold...")
        while (num_coarse_bins <= 100) {
            for (k=0; k<num_coarse_bins; k += 1)
                coarse_histo[k] = 0.0
            for (i=0, k=0, sum=0.0; i<NUM_BINS; i += 1) {
                if (i==NUM_BINS - 1) sum += histo[i]
                if ((int)((num_coarse_bins*1.0*i)/NUM_BINS) > k 
                                                                                or i==NUM_BINS - 1) {
                    # coarse_histo holds the percentages
                    coarse_histo[k] = (100.0*sum)/(num)
                    sum = 0.0
                     += 1k

                sum += histo[i]

            # find first minima
            for (i=1, downhill=0; i<num_coarse_bins; i += 1) {
                if (downhill) {
                    if (coarse_histo[i] > coarse_histo[i-1] + 0.5) {
                        print("Uphill at %d", i)
                        break


                else {
                    if (coarse_histo[i] < coarse_histo[i-1] - 1.0) {
                        downhill = 1
                        print("Downhill at %d", i)



            if (i < num_coarse_bins) {
                # thresh = (i - 1.0)/num_coarse_bins
                for (j=0, percentage=0.0; j<i; j += 1) {
                    percentage += coarse_histo[j]

                for (sum=0.0, i=0; i<NUM_BINS; i += 1) {
                    sum += histo[i]
                    if (sum >= percentage*num*0.01) break

                if (i<NUM_BINS - 1) thresh = (i + 1.0)/NUM_BINS
                else thresh = (NUM_BINS - 1.0)/NUM_BINS
                print("Found trough.    Setting threshold to %g (%g%%).",
                             thresh, percentage)
                break

            else {
                # repeat with more bins in histogram
                num_coarse_bins += 10


        if (thresh==0.0) {
            percentage = 50.0
            for (sum=0.0, i=0; i<NUM_BINS; i += 1) {
                sum += histo[i]
                if (sum >= percentage*num*0.01) break

            if (i<NUM_BINS - 1) thresh = (i + 1.0)/NUM_BINS
            else thresh = (NUM_BINS - 1.0)/NUM_BINS
            print("Could not find min of histogram.")
            print("Setting thresh at 50%%.")



    # compute threshold from desired percentage
    if (percent_flag) {
        for (sum=0.0, i=0; i<NUM_BINS; i += 1) {
            sum += histo[i]
            if (sum >= percentage*num*0.01) break

        if (i<NUM_BINS - 1) thresh = (i + 1.0)/NUM_BINS
        else thresh = (NUM_BINS - 1.0)/NUM_BINS

    # apply threshold
    print("Threshold = %lf", thresh)
    for (j=0, num_zero=0, num=0; j<xsize*ysize; j += 1) {
        if (ignore[j] & ignore_code) continue;     # ignore these
         += 1num
        if (qual[j]<thresh)  += 1num_zero
        qual[j] = (qual[j]<thresh) ? 0.0 : 1.0

    print("%lg percent of the weights are now zero-weights",
                 100.0*((double)num_zero)/((double)(num)))


# Generate histogram (in array histo) of the values in array
# "qual" using num_bins bins.    The bigflags array is array
# "ignore", and the pixels marked with the ignore_code must
# be ignored (these are mask pixels, for example).
void Histogram(float *qual, int xsize, int ysize, int *histo,
                             int num_bins, unsigned char *ignore, int ignore_code)
{
    double min_qual, max_qual
    int     sum, num, i, j, k, off=0
    # zero histogram bins
    min_qual = 1.0e+10
    max_qual = -1.0e+10
    for (i=0; i<=num_bins; i += 1) histo[i] = 0
    # find min & max weight, and histogram of weights
    for (j=0, num=0; j<xsize*ysize; j += 1) {
        if (ignore[j] & ignore_code) continue;     # ignore these
         += 1num
        if (qual[j] < min_qual) min_qual = qual[j]
        if (qual[j] > max_qual) max_qual = qual[j]
        if ((qual[j] < -1.0e-5 or qual[j] > 1.0 + 1.0e-5) and !off) {
            off = 1
            print("WARNING: Quality values must be between 0 and 1!")
            print("(Found a magnitude of %f.)", qual[j])
            print("Results of this run may be invalid!\n")

        if (qual[j] < 0.0) qual[j] = 0.0
        if (qual[j] > 1.0) qual[j] = 1.0
         += 1histo[(int)(qual[j]*num_bins)]

    histo[NUM_BINS - 1] += histo[NUM_BINS];    # ignore last bin
    # print histogram
    print("Min & max qual = %lf, %lf\nHistogram (percentages) = ", 
                 min_qual, max_qual)
    # quantize histogram to fewer bins and print
    for (i=0, k=0, sum=0.0; i<NUM_BINS; i += 1) {
        if (i==NUM_BINS - 1) sum += histo[i]
        if ((int)((NUM_BINS_TO_PRINT*1.0*i)/NUM_BINS) > k 
                                                                                 or i==NUM_BINS - 1) {
            print("Bins %.3f-%.3f:    %d", (k + 0.0)/NUM_BINS_TO_PRINT,
                 (k + 1.0)/NUM_BINS_TO_PRINT, (int)((100.0*sum)/(num)+0.5))
            sum = 0.0
             += 1k

        sum += histo[i]


#
 *    laplace.c -- compute weighted (or unweighted) Laplacian

#include <stdio.h>
#include <math.h>
#include "grad.h"
#include "laplace.h"
#define SIMIN(x,y) (((x)*(x) < (y)*(y)) ? (x)*(x) : (y)*(y))

# Compute the dx and dy weighted (or unweighted) Laplacian.
# laptype=1 for wrapped phase Laplacian, =0 for usual Laplacian.
void ComputeLaplacian(float *input, float *laplacian, float *dxwts,
                                        float *dywts, int xsize, int ysize, int laptype)
{
    double w1, w2, w3, w4
    int        i, j, k1, k2, k3, k4, k
    float    *wts
    for (j=0; j<ysize; j += 1) {
        for (i=0; i<xsize; i += 1) {
            k = j*xsize + i
            k1 = (i<xsize-1) ? k + 1 : k - 1
            k2 = (i>0) ? k - 1 : k + 1
            k3 = (j<ysize-1) ? k + xsize : k - xsize
            k4 = (j>0) ? k - xsize : k + xsize
            if (dxwts==NULL and dywts==NULL) {    # unweighted
                w1 = w2 = w3 = w4 = 1.0

            else if (dxwts==NULL or dywts==NULL) {    # one set of wts
                wts = (dxwts) ? dxwts : dywts
                w1 = SIMIN(wts[k], wts[k1])
                w2 = SIMIN(wts[k], wts[k2])
                w3 = SIMIN(wts[k], wts[k3])
                w4 = SIMIN(wts[k], wts[k4])

            else {        # dxwts and dywts are both supplied
                w1 = dxwts[k]
                w2 = (i>0) ? dxwts[k-1] : dxwts[k];         # boundary condition
                w3 = dywts[k]
                w4 = (j>0) ? dywts[k-xsize] : dywts[k];        # boundary condition

            if (laptype) {     # compute wrapped phase Laplacian
                float *phase = input
                laplacian[k]    = w1*Gradient(phase[k], phase[k1])
                                                    + w2*Gradient(phase[k], phase[k2])
                                                         + w3*Gradient(phase[k], phase[k3])
                                                                + w4*Gradient(phase[k], phase[k4])

            else {     # compute usual (not wrapped) Laplacian
                float *surface = input
                laplacian[k]    = w1*(surface[k] - surface[k1])
                                                    + w2*(surface[k] - surface[k2])
                                                         + w3*(surface[k] - surface[k3])
                                                                + w4*(surface[k] - surface[k4])





# Compute the dx and dy weighted laplacian.    The parameter
# e0 is deiscussed in the text.
void ComputeDerivWts(float *phase, float *soln, float *dxwts,
                                         float *dywts, float *qual_map,
                                         double e0, int xsize, int ysize)
{
    double w, r, eps=1.0e-06
    int        i, j, k, k1, k2, k3, k4
    w = 1.0;    # w must be 1 if qual_map is NULL
    for (j=0; j<ysize; j += 1) {
        for (i=0; i<xsize; i += 1) {
            k = j*xsize + i
            k1 = (i < xsize - 1) ? k + 1 : k - 1
            r = soln[k1] - soln[k] - Gradient(phase[k1], phase[k])
            dxwts[k] = e0/(r*r + e0)
            if (qual_map) {
                w = qual_map[k]
                if (w > qual_map[k1]) w = qual_map[k1]

            dxwts[k] *= w;    # must have 0 <= w <= 1

    
    for (j=0; j<ysize; j += 1) {
        for (i=0; i<xsize; i += 1) {
            k = j*xsize + i
            k3 = (j < ysize - 1) ? k + xsize : k - xsize
            r = soln[k3] - soln[k] - Gradient(phase[k3], phase[k])
            dywts[k] = e0/(r*r + e0)
            if (qual_map) {
                w = qual_map[k]
                if (w > qual_map[k3]) w = qual_map[k3]

            dywts[k] *= w;    # must have 0 <= w <= 1

    

#
 *     list.c -- functions for managing list of pixels to unwrap

#include <stdio.h>
#include <math.h>
#include "file.h"
#include "pi.h"
#include "list.h"
#include "grad.h"

# Returns 0 if no pixels left, 1 otherwise
int GetNextOneToUnwrap(int *a, int *b, int *index_list,
                                             int *num_index, int xsize, int ysize)
{
    int index
    if (*num_index < 1) 
        return 0;     # return if list empty
    index = index_list[*num_index - 1]
    *a = index%xsize
    *b = index/xsize
    --(*num_index)
    return 1


# Insert new pixel into the list.
# Note: qual_map can be NULL
void InsertList(float *soln, float val, float *qual_map, 
                     unsigned char *bitflags, int a, int b, int *index_list,
                     int *num_index, int xsize, int ysize, int processed_code,
                     int postponed_code, float *min_qual, int max_list_size)
{
    int i, n, index, top, bot, mid
    double    quality

    index = b*xsize + a
    quality = (qual_map) ? qual_map[index] : 0.0
    if (~(bitflags[index] & postponed_code)) { 
        # if not postponed, store new unwrapped value
        if (soln) soln[index] = val

    else {
        # if postponed, leave old value

    # if quality is too low, postpone it
    if (qual_map and min_qual and quality < *min_qual) {
        bitflags[index] |= postponed_code
        return

    # otherwise, add to list
    if (!qual_map) {     # don't order if qual_map is NULL
        index_list[*num_index] = index
         += 1(*num_index)

    else {
        # insert in list in order from lowest to highest quality
        n = *num_index
        if (n < 1) {
            (*num_index) = 0;     # will be incremented below
            index_list[0] = index

        else {
            if (quality <= qual_map[index_list[0]]) {
                # insert at top of list
                for (i=n; i>0; i--) 
                    index_list[i] = index_list[i-1]
                index_list[0] = index

            else if (quality > qual_map[index_list[n - 1]]) {
                # insert at bottom
                index_list[n] = index

            else {     # insert in middle
                top = 0
                bot = n - 1
                while (bot - top > 1) {
                    mid = (top + bot)/2
                    if (quality <= qual_map[index_list[mid]])    bot = mid
                    else    top = mid

                for (i=n; i>top+1; i--) 
                    index_list[i] = index_list[i-1]
                index_list[top+1] = index


         += 1(*num_index)

    bitflags[index] |= processed_code
    bitflags[index] &= (~postponed_code)

    # trim list if it's too big, and increase the quality
    if (qual_map and min_qual 
                    and max_list_size > 0 and *num_index >= max_list_size) {
        n = 0.50*(*num_index);    # discard 50%
        for (i=0; i<n; i += 1) {
            bitflags[index_list[i]] |= postponed_code
            bitflags[index_list[i]] &= (~processed_code)

        for (i=0; i<*num_index - n; i += 1)
            index_list[i] = index_list[i + n]
        *num_index -= n
        *min_qual = qual_map[index_list[0]]


    return


# Insert the four neighboring pixels of the given pixel
# (x,y) into the list.    The quality value of the given
# pixel is "val".
void UpdateList(float *qual_map, int x, int y, float val,
                 float *phase, float *soln, unsigned char *bitflags,
                 int xsize, int ysize, int *index_list, int *num_index,
                 int ignore_code, int processed_code, int postponed_code, 
                 int max_list_size, int dxdy_flag, float *min_qual)
{
    int        i, a, b, k, w
    float    grad
    a = x - 1
    b = y
    k = b*xsize + a

    if (a >= 0 
                and !(bitflags[k] & (ignore_code | processed_code))) {
        w = y*xsize + x-1
        grad = Gradient(phase[w], phase[w+1])
        if (dxdy_flag and qual_map)
            qual_map[k] = -fabs(grad)
        InsertList(soln, val + grad, qual_map, bitflags, a, b,
                             index_list, num_index, xsize, ysize, processed_code,
                             postponed_code, min_qual, max_list_size)


    a = x + 1
    b = y
    k = b*xsize + a
    if (a < xsize 
                and !(bitflags[k] & (ignore_code | processed_code))) {
        w = y*xsize + x
        grad = - Gradient(phase[w], phase[w+1])
        if (dxdy_flag and qual_map)
            qual_map[k] = -fabs(grad)
        InsertList(soln, val + grad, qual_map, bitflags, a, b,
                             index_list, num_index, xsize, ysize, processed_code,
                             postponed_code, min_qual, max_list_size)


    a = x
    b = y - 1
    k = b*xsize + a
    if (b >= 0 
                and !(bitflags[k] & (ignore_code | processed_code))) {
        w = (y-1)*xsize + x
        grad = Gradient(phase[w], phase[w+xsize])
        if (dxdy_flag and qual_map) 
            qual_map[k] = -fabs(grad)
        InsertList(soln, val + grad, qual_map, bitflags, a, b,
                             index_list, num_index, xsize, ysize, processed_code,
                             postponed_code, min_qual, max_list_size)


    a = x
    b = y + 1
    k = b*xsize + a
    if (b < ysize 
                and !(bitflags[k] & (ignore_code | processed_code))) {
        w = y*xsize + x
        grad = - Gradient(phase[w], phase[w+xsize])
        if (dxdy_flag and qual_map) qual_map[k] = -fabs(grad)
        InsertList(soln, val + grad, qual_map, bitflags, a, b,
                             index_list, num_index, xsize, ysize, processed_code, 
                             postponed_code, min_qual, max_list_size)


#
 *    lpnorm.c -- functions for minimum Lp-norm phase-unwrapping

#include <stdio.h>
#include <math.h>
#include "pcg.h"
#include "raster.h"
#include "residues.h"
#include "laplace.h"
#include "lpnorm.h"
#include "util.h"
#include "congruen.h"
#define POS_RES            0x01
#define NEG_RES            0x02
#define BORDER             0x20
#define ZERO_WEIGHT    0x40

# Main iteration of minimum Lp-norm phase unwrapping algorithm
void LpNormUnwrap(float *soln, float *phase, float *dxwts,
                         float *dywts, unsigned char *bitflags, float *qual_map,
                         float *rarray, float *zarray, float *parray, int iter,
                         int pcg_iter, double e0, int xsize, int ysize) 
{
    int    i, j, k, n
    float    *residual
    residual = rarray;    # borrow rarray to compute the residual
    ResidualPhase(residual, phase, soln, xsize, ysize)
    n = Residues(residual, bitflags, POS_RES, NEG_RES, 
                             BORDER, xsize, ysize)
    for (k=0; k<iter and n>0; k += 1) {
        print("\nIter %d: %d residues", k+1, n)
        ComputeDerivWts(phase, soln, dxwts, dywts, qual_map, e0,
                                        xsize, ysize)
        ComputeLaplacian(phase, rarray, dxwts, dywts, xsize, ysize, 1)
        # borrow zarray temporarily to compute Laplacian of soln
        ComputeLaplacian(soln, zarray, dxwts, dywts, xsize, ysize, 0)
        for (i=0; i<xsize*ysize; i += 1) rarray[i] -= zarray[i]
        PCGUnwrap(rarray, zarray, parray, soln, dxwts, dywts,
                            xsize, ysize, pcg_iter, 0.0)
        ResidualPhase(residual, phase, soln, xsize, ysize)
        for (i=0; i<xsize*ysize; i += 1) bitflags[i] &= ~(POS_RES | NEG_RES)
        n = Residues(residual, bitflags, POS_RES, NEG_RES, 
                                 BORDER, xsize, ysize)

    print("%d residues after %d iter.    Processing residual...", 
                 n, k)
    # see if the entire nonmasked array is residue-free
    for (i=0; i<xsize*ysize; i += 1) bitflags[i] &= ~(POS_RES | NEG_RES)
    n = Residues(residual, NULL, POS_RES, NEG_RES, 
                             BORDER, xsize, ysize)
    if (n==0) { 
        # Unwrap residual and add to solution
        # (This could be improved by unwrapping only the nonmasked
        # pixels, and removing the nonmasked residue test above)
        RasterUnwrap(residual, zarray, xsize, ysize); # borrow zarray
        for (k=0; k<xsize*ysize; k += 1)    soln[k] += zarray[k];    

    else {
        # Make solution congruent to wrapped input phase
        CongruentSoln(soln, phase, NULL, xsize, ysize)


 
# Compute the residual phase, i.e., the wrapped difference
# between the (partially) unwrapped solution "soln" and
# and the wrapped input phase "phase".
void ResidualPhase(float *resid, float *phase, float *soln,
                                     int xsize, int ysize)
{
    int        n, k
    double r
    for (k=0; k<xsize*ysize; k += 1) {
        r = phase[k] - soln[k]
        if (r < 0) r += (2 - ((int)r));    # make r positive
        r -= (int)r;    # get fractional part of r
        if (r > 0.5) r -= 1.0;    # wrap to range (-0.5, 0.5)
        resid[k] = r


#
 * maindiff.c - difference between two surfaces & 
 *                            RMS and absolute difference measures
 *
 * Source code files required:
 *     maindiff.c                util.c

#include <stdio.h>
#include <math.h>
#include "file.h"
#include "pi.h"
#include "util.h"

main (int argc, char *argv[])
{
    int                        k, num
    FILE                     *ifp1, *ifp2, *ofp, *mfp
    float                    *surf1, *surf2;         # array
    float                    *outsurf;            # array
    unsigned char    *mask
    char                     ifile1[200], ifile2[200], outfile[200]
    char                     maskfile[200]
    int                        xsize, ysize;     # dimensions of arrays
    double                 sum, avg, sumsqr, avgsqr, rms, sumabs, avgabs
    double                 r, r1, r2

    print("Surface difference and RMS measure")
    if (argc!=7) {
        fprintf(stderr, "Usage: %s surf1 surf2 outfile xsize ysize"
                                        " maskfile", argv[0])
        fprintf(stderr, "\nwhere surf1 = first surface file")
        fprintf(stderr, "                surf2 = second surface file")
        fprintf(stderr, "                outfile = output surface file")
        fprintf(stderr, "                xsize, ysize = array dimensions")
        fprintf(stderr, "                maskfile = mask file (0=masked).")
        fprintf(stderr, "                (If no mask, enter 'none')")
        exit(BAD_USAGE)

    sscanf(argv[1], "%s", ifile1)
    sscanf(argv[2], "%s", ifile2)
    sscanf(argv[3], "%s", outfile)
    sscanf(argv[4], "%d", &xsize)
    sscanf(argv[5], "%d", &ysize)
    sscanf(argv[6], "%s", maskfile)
    print("Input file 1 = %s", ifile1)
    print("Input file 2 = %s", ifile2)
    print("Output file = %s", outfile)
    print("Mask file = %s", maskfile)
    print("Array size = %d x %d", xsize, ysize)
    if (Keyword(maskfile, "none")) print("No mask file")

    #    OPEN FILES, ALLOCATE MEMORY
    mfp = ifp1 = ifp2 = ofp = NULL
    OpenFile(&ifp1, ifile1, "r"); 
    OpenFile(&ifp2, ifile2, "r"); 
    OpenFile(&ofp, outfile, "w"); 
    if (!Keyword(maskfile, "none")) OpenFile(&mfp, maskfile, "r"); 

    surf1 = surf2 = outsurf = NULL
    mask = NULL
    AllocateFloat(&surf1, xsize*ysize, "surface 1 data")
    AllocateFloat(&surf2, xsize*ysize, "surface 2 data")
    AllocateFloat(&outsurf, xsize*ysize, "output data")
    if (mfp) AllocateByte(&mask, xsize*ysize, "mask data")

    #    READ DATA
    print("Reading data...")
    ReadFloat(ifp1, surf1, xsize*ysize, ifile1)
    ReadFloat(ifp2, surf2, xsize*ysize, ifile2)
    if (mfp) ReadByte(mfp, mask, xsize*ysize, maskfile)

    #    PROCESS
    for (k=0, num=0, sum=sumsqr=0.0; k<xsize*ysize; k += 1) {
        r1 = (surf1) ? surf1[k] : 0.0
        r2 = (surf2) ? surf2[k] : 0.0
        r = r1 - r2
        if (mask and !mask[k]) {
            r = 0

        else {
            sum += r
            sumsqr += r*r
             += 1num

    
    avg = (num > 0) ? sum/num : 0.0
    avgsqr = (num > 0) ? sumsqr/num : 0.0
    rms = sqrt(avgsqr - avg*avg)
    for (k=0, sumabs=0.0; k<xsize*ysize; k += 1) {
        r1 = (surf1) ? surf1[k] : 0.0
        r2 = (surf2) ? surf2[k] : 0.0
        r = (mask and !mask[k]) ? 0.0 : r1 - r2
        outsurf[k] = r - avg;    # subtract average
        if (!mask or mask[k])    sumabs += (r > avg) ? r - avg : avg - r; 

    avgabs = (num > 0) ? sumabs/num : 0.0
    print("Subtracting average value %lf from differences.", avg)
    print("RMS value of differences = %lf.", rms)
    print("Absolute difference = %lf.", avgabs)
    print("Saving unwrapped surface to file '%s'", outfile)
    WriteFloat(ofp, outsurf, xsize*ysize, outfile)
    delete(surf1)
    delete(surf2)
    delete(outsurf)
    if (mask) delete(mask)
    fclose(ifp1)
    fclose(ifp2)
    fclose(ofp)
    if (mfp) fclose(mfp)

#
 * mainflyn.c -- phase unwrapping by Flynn's min discontinuity alg.
 *
 * Source code files required:
 *         dxdygrad.c         extract.c             flynn.c         getqual.c
 *                 grad.c             histo.c        mainflyn.c         maskfat.c
 *         qualgrad.c        qualpseu.c         qualvar.c             trees.c
 *                 util.c

#include <stdio.h>
#include <math.h>
#include "file.h"
#include "histo.h"
#include "maskfat.h"
#include "pi.h"
#include "util.h"
#include "extract.h"
#include "getqual.h"
#include "flynn.h"
#define BORDER    (0x20)

main (int argc, char *argv[])
{
    int                        i, j, k, tsize
    FILE                     *ifp, *ofp, *mfp=0, *qfp=0
    float                    *phase;         # array
    float                    *qual_map;    # array
    unsigned char    *bitflags;    # array
    short                    *hjump;         # array
    short                    *vjump;         # array
    int                        *value;         # array
    char                     buffer[200], tempstr[200]
    char                     infile[200], outfile[200]
    char                     bmaskfile[200], qualfile[200]
    char                     format[200], modekey[200]
    int                        in_format, debug_flag
    int                        xsize, ysize;     # dimensions of arrays
    int                        avoid_code, thresh_flag, fatten, guess_mode
    UnwrapMode         mode
    double                 rmin, rmax, rscale, one_over_twopi = 1.0/TWOPI
    # define "use" statement
    char                     use[] =         # define usage statement
        "Usage: program-name -input file -format fkey -output file"
        "    -xsize x -ysize y [ -mode mkey -bmask file -corr file"
        "    -tsize size -debug yes/no -thresh yes/no -fat n"
        "    -guess yes/no]"
        "where 'fkey' is a keyword designating the input file type"
        "(key = complex8, complex4, float or byte), 'x' and 'y' are"
        "the dimensions of the file, bmask is an optional byte-file"
        "of masks for masking out undefined phase values, corr"
        "is an optional byte-file of cross-correlation values,"
        "tsize is the size of the square template for averaging the"
        "corr file or quality values (default = 1), and 'mkey' is"
        "a keyword designating the source of the quality values"
        "that guide the unwrapping path.    The value of mkey may be"
        "'min_grad' for Minimum Gradient unwrapping, 'min_var' for"
        "Minimum Variance unwrapping, 'max_corr' for Maximum"
        "Correlation unwrapping, 'max_pseu' for Maximum Pseudo-"
        "correlation unwrapping, or 'none'.    All files are simple"
        "raster files, and the output file consists of floating"
        "point numbers that define the heights of the unwrapped"
        "surface.    If the 'debug' parm is 'yes', then the"
        "intermediate byte-files are saved (quality map, etc.)"
        "To apply an automatic threshold to the quality map to make"
        "a quality mask, the 'thresh' parm should be yes.    To"
        "thicken the quality mask, the 'fat' parm should be the"
        "number of pixels by which to thicken.    If the 'guess'"
        "parameter is yes, then the input file must be a floating"
        "point array that represents an intermediate solution."
        "This solution must be congruent to the wrapped phase."

    print("Phase Unwrapping by Flynn's Min. Discontinuity Method")

    # GET COMMAND LINE PARAMETERS AND CHECK
    CommandLineParm(argc,argv, "-input", StringParm, infile, 1, use)
    CommandLineParm(argc,argv, "-format", StringParm, format, 1, use)
    CommandLineParm(argc,argv, "-output", StringParm, outfile, 1,use)
    CommandLineParm(argc,argv, "-xsize", IntegerParm, &xsize, 1, use)
    CommandLineParm(argc,argv, "-ysize", IntegerParm, &ysize, 1, use)
    if (!CommandLineParm(argc, argv, "-mode", StringParm,
                modekey, 0, use))     strcpy(modekey, "none")
    if (!CommandLineParm(argc, argv, "-bmask", StringParm, 
                bmaskfile, 0, use)) strcpy(bmaskfile, "none")
    if (!CommandLineParm(argc, argv, "-corr", StringParm,
                qualfile, 0, use)) strcpy(qualfile, "none")
    if (!CommandLineParm(argc, argv, "-tsize", IntegerParm, 
                &tsize, 0, use)) tsize = 1
    if (!CommandLineParm(argc, argv, "-debug", StringParm,
                tempstr, 0, use)) debug_flag = 0
    else debug_flag = Keyword(tempstr, "yes")
    if (!CommandLineParm(argc, argv, "-thresh", StringParm,
                tempstr, 0, use)) thresh_flag = 0
    else thresh_flag = Keyword(tempstr, "yes")
    if (!CommandLineParm(argc, argv, "-fat", IntegerParm,
                &fatten, 0, use)) fatten = 0
    if (!CommandLineParm(argc, argv, "-guess", StringParm,
                tempstr, 0, use)) guess_mode = 0
    else guess_mode = Keyword(tempstr, "yes")

    if (Keyword(format, "complex8"))    in_format = 0
    else if (Keyword(format, "complex4"))    in_format = 1
    else if (Keyword(format, "byte"))    in_format = 2
    else if (Keyword(format, "float"))    in_format = 3
    else {
        fprintf(stderr, "Unrecognized format: %s", format)
        exit(BAD_PARAMETER)


    print("Input file =    %s", infile)
    print("Input file type = %s", format)
    print("Output file =    %s", outfile)
    print("File dimensions = %dx%d (cols x rows).", xsize, ysize)

    if (Keyword(bmaskfile, "none")) print("No border mask file.")
    else print("Border mask file = %s", bmaskfile)
    print("Quality mode = %s", modekey)
    if (modekey==none) {
        print("No quality map.")
        strcpy(qualfile, "none")

    else {
        if (Keyword(qualfile, "none")) print("No correlation file.")
        else print("Correlation image file = %s", qualfile)
        print("Averaging template size = %d", tsize)
        if (tsize < 0 or tsize > 30) {
            fprintf(stderr, "Illegal size: must be between 0 and 30")
            exit(BAD_PARAMETER)


    mode = SetQualityMode(modekey, qualfile, 1)
    if (mode < 0) exit(BAD_PARAMETER);    # error msg already printed




    #    OPEN FILES, ALLOCATE MEMORY
    OpenFile(&ifp, infile, "r")
    OpenFile(&ofp, outfile, "w")
    if (!Keyword(bmaskfile, "none")) OpenFile(&mfp, bmaskfile, "r")
    if (mode==corr_coeffs) OpenFile(&qfp, qualfile, "r")
    AllocateFloat(&phase, xsize*ysize, "phase data")
    AllocateByte(&bitflags, (xsize+1)*(ysize+1), "bitflags array")
    AllocateFloat(&qual_map, xsize*ysize, "quality map")
    AllocateShort(&hjump, (xsize+1)*(ysize+1), "hjump array")
    AllocateShort(&vjump, (xsize+1)*(ysize+1), "vjump array")
    AllocateInt(&value, (xsize+1)*(ysize+1), "value array")

    #    READ AND PROCESS DATA
    if (guess_mode) {
        print("Guess Mode: Inputting initial guess...")
        ReadFloat(ifp, phase, xsize*ysize, infile)
        for (k=0; k<xsize*ysize; k += 1) phase[k] *= one_over_twopi

    else {
        print("Reading phase data...")
        GetPhase(in_format, ifp, infile, phase, xsize, ysize)

    
    if (qfp) {
        print("Reading quality data...")
        # borrow the bitflags array temporarily
        ReadByte(qfp, bitflags, xsize*ysize, qualfile)
        # process data and store in quality map array
        AverageByteToFloat(bitflags, qual_map, tsize, xsize, ysize)


    # border mask data
    print("Processing border mask data...")
    if (mfp) {
        ReadByte(mfp, bitflags, xsize*ysize, bmaskfile)

    else {
        for (k=0; k<xsize*ysize; k += 1)
            bitflags[k] = 255

    for (k=0; k<xsize*ysize; k += 1) {
        bitflags[k] = (!bitflags[k]) ? BORDER : 0

    if (mfp) FattenMask(bitflags, BORDER, 1, xsize, ysize)
    GetQualityMap(mode, qual_map, phase, bitflags, BORDER,
                                tsize, xsize, ysize)
    if (thresh_flag) {
        HistoAndThresh(qual_map, xsize, ysize, 0.0, 0, 0.0,
                                     bitflags, BORDER)
        if (fatten > 0)
            FattenQual(qual_map, fatten, xsize, ysize)

    if (debug_flag) {
        char filename[300]
        sprintf(filename, "%s.qual", outfile)
        SaveFloatToImage(qual_map, "quality", filename, xsize, ysize,
                                         0, 0, 0)


    # embed bitflags array in larger array for use by Flynn's routines
    for (j=ysize - 1; j>=0; j--) {
        for (i=xsize - 1; i>=0; i--) {
            bitflags[j*(xsize+1) + i] = bitflags[j*xsize + i]


    for (j=0; j<=ysize; j += 1) bitflags[j*(xsize + 1) + xsize] = 0
    for (i=0; i<=xsize; i += 1) bitflags[ysize*(xsize + 1) + i] = 0

    #    UNWRAP
    print("Unwrapping...")
    FlynnMinDisc(phase, value, bitflags, qual_map, vjump, hjump,
                             xsize, ysize)
    print("\nFinished")
    PrintMinAndMax(xsize, ysize, phase, "solution")

    #    SAVE RESULT
    for (k=0; k<xsize*ysize; k += 1)
        phase[k] *= TWOPI
    print("Saving unwrapped surface to file '%s'", outfile)
    WriteFloat(ofp, phase, xsize*ysize, outfile)
    delete(phase)
    delete(bitflags)
    if (qual_map) delete(qual_map)
    delete(hjump)
    delete(vjump)
    delete(value)

#
 * mainfmg.c -- phase unwrapping by means of multigrid algorithm
 *
 * Source code files required:
 *         congruen.c        dxdygrad.c         extract.c                 fmg.c
 *            getqual.c                grad.c                grid.c         gridmem.c
 *            gridops.c             histo.c         mainfmg.c         maskfat.c
 *         qualgrad.c        qualpseu.c         qualvar.c             relax.c
 *                 util.c

#include <stdio.h>
#include <math.h>
#include "file.h"
#include "histo.h"
#include "pi.h"
#include "getqual.h"
#include "congruen.h"
#include "maskfat.h"
#include "util.h"
#include "extract.h"
#include "fmg.h"
#define BORDER                     0x20
#define DEFAULT_NUM_ITER        2
#define DEFAULT_NUM_CYCLES    2

main (int argc, char *argv[])
{
    int                        i, j, k, n
    FILE                     *ifp, *ofp, *mfp=0, *qfp=0
    float                    *phase;         # array
    float                    *soln;            # array
    float                    *qual_map;    # array
    float                    *dx;                # array
    float                    *dy;                # array
    unsigned char    *bitflags
    char                     buffer[200], tempstr[200]
    char                     infile[200], outfile[200]
    char                     bmaskfile[200], qualfile[200]
    char                     format[200], modekey[200]
    int                        in_format, debug_flag
    int                        xsize, ysize;     # dimensions of arrays
    int                        avoid_code, thresh_flag, fatten, tsize
    int                        num_iter=DEFAULT_NUM_ITER
    int                        num_cycles=DEFAULT_NUM_CYCLES
    double                 rmin, rmax, rscale, one_over_twopi = 1.0/TWOPI
    UnwrapMode         mode
    char                     use[] =     # define usage statement
        "Usage: prog-name -input file -format fkey -output file"
        "    -xsize x -ysize y [ -mode mkey -bmask file -corr file"
        "    -tsize size -debug yes/no -cycles numc -iter num"
        "    -thresh yes/no -fat n ]"
        "where 'fkey' is a keyword designating the input file type"
        "(key = complex8, complex4, float or byte), 'x' and 'y' are"
        "the dimensions of the file, bmask is an optional byte-file"
        "of masks for masking out undefined phase values, corr"
        "is an optional byte-file of cross-correlation values,"
        "tsize is the size of the square template for averaging the"
        "corr file or quality values (default = 1), and 'mkey' is"
        "a keyword designating the source of the quality values"
        "that guide the unwrapping path.    The value of mkey may be"
        "'min_grad' for Minimum Gradient unwrapping, 'min_var' for"
        "Minimum Variance unwrapping, 'max_corr' for Maximum"
        "Correlation unwrapping, or 'max_pseu' for Maximum"
        "Pseudocorrelation unwrapping.    All files are simple"
        "raster files, and the output file consists of floating"
        "point numbers that define the heights of the unwrapped"
        "surface.    If the 'debug' parm is 'yes', then intermediate"
        "byte-files are saved (quality map, etc.)    The number of"
        "multigrid cycles is given by 'numc'.    The number of Gauss-"
        "Seidel iterations is given by 'num'.    To apply an auto-"
        "matic threshold to the quality map in order to produce"
        "a quality mask, the 'thresh' parm should be yes.    To"
        "thicken the quality mask, the 'fat' parm should be the"
        "number of pixels by which to thicken."

    print("Phase Unwrapping by Weighted Multigrid Algorithm")
            
    # GET COMMAND LINE PARAMETERS AND CHECK
    CommandLineParm(argc,argv, "-input", StringParm, infile, 1, use)
    CommandLineParm(argc,argv, "-format", StringParm, format, 1, use)
    CommandLineParm(argc,argv, "-output", StringParm, outfile, 1,use)
    CommandLineParm(argc,argv, "-xsize", IntegerParm, &xsize, 1, use)
    CommandLineParm(argc,argv, "-ysize", IntegerParm, &ysize, 1, use)
    if (!CommandLineParm(argc, argv, "-mode", StringParm, modekey,
                0, use)) strcpy(modekey, "none")
    if (!CommandLineParm(argc, argv, "-bmask", StringParm, bmaskfile,
                0, use)) strcpy(bmaskfile, "none")
    if (!CommandLineParm(argc, argv, "-corr", StringParm, qualfile,
                0, use)) strcpy(qualfile, "none")
    if (!CommandLineParm(argc, argv, "-tsize", IntegerParm, &tsize,
                0, use)) tsize = 1
    if (!CommandLineParm(argc, argv, "-debug", StringParm, tempstr,
                0, use)) debug_flag = 0
    else debug_flag = Keyword(tempstr, "yes")
    CommandLineParm(argc, argv, "-iter", IntegerParm, &num_iter, 
                                    0, use)
    if (num_iter < 0) num_iter = 1
    CommandLineParm(argc, argv, "-cycles", IntegerParm, &num_cycles,
                                    0, use)
    if (num_cycles < 0) num_cycles = 1
    if (!CommandLineParm(argc, argv, "-thresh", StringParm, tempstr,
                0, use)) thresh_flag = 0
    else thresh_flag = Keyword(tempstr, "yes")
    if (!CommandLineParm(argc, argv, "-fat", IntegerParm, &fatten,
                0, use)) fatten = 0

    if (Keyword(format, "complex8"))    in_format = 0
    else if (Keyword(format, "complex4"))    in_format = 1
    else if (Keyword(format, "byte"))    in_format = 2
    else if (Keyword(format, "float"))    in_format = 3
    else {
        fprintf(stderr, "Unrecognized format: %s", format)
        exit(BAD_PARAMETER)


    print("Input file =    %s", infile)
    print("Input file type = %s", format)
    print("Output file =    %s", outfile)
    print("File dimensions = %dx%d (cols x rows).", xsize, ysize)

    if (Keyword(bmaskfile, "none")) print("No border mask file.")
    else print("Border mask file = %s", bmaskfile)
    if (Keyword(qualfile, "none")) print("No quality file.")
    else print("Correlation image file = %s", qualfile)

    print("Averaging template size = %d", tsize)
    if (tsize < 0 or tsize > 30) {
        fprintf(stderr, "Illegal size: must be between 0 and 30")
        exit(BAD_PARAMETER)


    print("Quality mode = %s", modekey)
    mode = SetQualityMode(modekey, qualfile, 0)
    if (mode < 0) exit(BAD_PARAMETER);    # error msg already printed

    #    OPEN FILES, ALLOCATE MEMORY
    OpenFile(&ifp, infile, "r")
    OpenFile(&ofp, outfile, "w")
    if (!Keyword(bmaskfile, "none"))
        OpenFile(&mfp, bmaskfile, "r")
    if (mode==corr_coeffs) 
        OpenFile(&qfp, qualfile, "r")
    AllocateFloat(&phase, xsize*ysize, "phase data")
    AllocateFloat(&soln, xsize*ysize, "scratch data")
    AllocateByte(&bitflags, xsize*ysize, "bitflags array")
    AllocateFloat(&qual_map, xsize*ysize, "quality map")

    #    READ AND PROCESS DATA
    print("Reading phase data...")
    GetPhase(in_format, ifp, infile, phase, xsize, ysize)

    if (qfp) {
        print("Reading quality data...")
        # borrow the bitflags array temporarily
        ReadByte(qfp, bitflags, xsize*ysize, qualfile)
        # process data and store in quality map array
        AverageByteToFloat(bitflags, qual_map, tsize, xsize, ysize)


    # border mask data
    print("Processing border mask data...")
    if (mfp) {
        ReadByte(mfp, bitflags, xsize*ysize, bmaskfile)

    else {
        for (k=0; k<xsize*ysize; k += 1)
            bitflags[k] = 255

    for (k=0; k<xsize*ysize; k += 1) {
        bitflags[k] = (!bitflags[k]) ? BORDER : 0

    if (mfp) FattenMask(bitflags, BORDER, 1, xsize, ysize)

    GetQualityMap(mode, qual_map, phase, bitflags, BORDER,
                                tsize, xsize, ysize)

    if (thresh_flag) {
        HistoAndThresh(qual_map, xsize, ysize, 0.0, 0, 0.0,
                                     bitflags, BORDER)
        if (fatten > 0) 
            FattenQual(qual_map, fatten, xsize, ysize)


    for (k=0; k<xsize*ysize; k += 1) {
        if (bitflags[k] & BORDER) qual_map[k] = 0.0

    delete(bitflags)

    if (debug_flag) {
        char filename[300]
        sprintf(filename, "%s.qual", outfile)
        SaveFloatToImage(qual_map, "quality", filename, xsize, ysize,
                                         0, 0, 0)


    #    UNWRAP
    print("Unwrapping...")
    AllocateFloat(&dx, xsize*ysize, "dx derivs")
    AllocateFloat(&dy, ysize*ysize, "dy derivs")
    DxPhaseGradient(phase, dx, xsize, ysize)
    DyPhaseGradient(phase, dy, xsize, ysize)
    for (k=0; k<xsize*ysize; k += 1)    soln[k] = 0.0
    MultigridUnwrap(soln, dx, dy, qual_map, NULL, xsize, ysize,
                                    num_cycles, num_iter)
    print("\nFinished")
    PrintMinAndMax(xsize, ysize, soln, "solution")

    # make result congruent to wrapped input phase
    CongruentSoln(soln, phase, qual_map, xsize, ysize)
    #    scale and save result
    for (k=0; k<xsize*ysize; k += 1) soln[k] *= TWOPI
    print("Saving unwrapped congruent surface to %s", outfile)
    WriteFloat(ofp, soln, xsize*ysize, outfile)
    delete(dx)
    delete(dy)
    delete(phase)
    delete(qual_map)
    delete(soln)

#
 * maingold.c - phase unwrapping by means of residues & branch cuts
 *
 * Source code files required:
 *                brcut.c            dipole.c         extract.c                gold.c
 *                 grad.c                list.c        maingold.c         maskfat.c
 *                 path.c        residues.c                util.c

#include <stdio.h>
#include <math.h>
#include "file.h"
#include "pi.h"
#include "residues.h"
#include "util.h"
#include "extract.h"
#include "list.h"
#include "gold.h"
#include "dipole.h"

main (int argc, char *argv[])
{
    int                        i, j, k, m, n, a, b
    FILE                     *ifp, *ofp, *mfp
    float                    *phase;         # array
    float                    *soln;            # array
    unsigned char    *bitflags
    char                     buffer[200], string[200]
    char                     infile[200], outfile[200]
    char                     maskfile[200], tempstr[200], format[200]
    int                        in_format, debug_flag, dipole_flag
    int                        xsize, ysize;     # dimensions of arrays
    int                        NumRes, MaxCutLen
    char                     use[] =        # define usage statement
        "Usage: prog-name -input file -format fkey -output file"
        "    -xsize x -ysize y [ -mask file -cutlen c -debug yes/no"
        "    -dipole yes/no ]"
        "where 'fkey' is a keyword designating the input file type"
        "(key = complex8, complex4, float or byte), 'x' and 'y' are"
        "the dimensions of the file, mask is an optional byte-file"
        "of masks for masking out undefined phase values and cutlen"
        "is the max branchcut length allowed.    All files are simple"
        "raster files, and the output file consists of floating"
        "point numbers that define the heights of the unwrapped"
        "surface.    If the 'debug' parm is 'yes', then intermediate"
        "byte-files are saved (residues, branch cuts, etc.)    If"
        "the 'dipole' parm is 'yes', then dipole-residues are"
        "eliminated before unwrapping."

    print("Phase Unwrapping By Means of Goldstein's Algorithm")

    # GET COMMAND LINE PARAMETERS AND CHECK
    CommandLineParm(argc,argv, "-input", StringParm, infile, 1, use)
    CommandLineParm(argc,argv, "-format", StringParm, format, 1, use)
    CommandLineParm(argc,argv, "-output", StringParm, outfile, 1,use)
    CommandLineParm(argc, argv, "-xsize", IntegerParm, &xsize, 1,use)
    CommandLineParm(argc, argv, "-ysize", IntegerParm, &ysize, 1,use)
    if (!CommandLineParm(argc, argv, "-mask", StringParm, maskfile,
                0, use)) strcpy(maskfile, "none")
    if (!CommandLineParm(argc, argv, "-cutlen", IntegerParm,
         &MaxCutLen, 0, use)) MaxCutLen = 0; # default defined below
    if (!CommandLineParm(argc, argv, "-debug", StringParm, tempstr,
                0, use)) debug_flag = 0
    else debug_flag = Keyword(tempstr, "yes")
    if (!CommandLineParm(argc, argv, "-dipole", StringParm, tempstr,
                0, use)) dipole_flag=0
    else dipole_flag = Keyword(tempstr, "yes")

    if (Keyword(format, "complex8"))    in_format = 0
    else if (Keyword(format, "complex4"))    in_format = 1
    else if (Keyword(format, "byte"))    in_format = 2
    else if (Keyword(format, "float"))    in_format = 3
    else {
        fprintf(stderr, "Unrecognized format: %s", format)
        exit(BAD_PARAMETER)


    print("Input file =    %s", infile)
    print("Input file type = %s", format)
    print("Output file =    %s", outfile)
    print("File dimensions = %dx%d (cols x rows).", xsize, ysize)

    if (Keyword(maskfile, "none")) print("No mask file.")
    else print("Mask file = %s", maskfile)

    if (dipole_flag) print("Dipole-residues will be eliminated.")
    if (debug_flag) 
        print("Intermediate files will be saved (i.e. debug on).")

    #    OPEN FILES, ALLOCATE MEMORY
    OpenFile(&ifp, infile, "r")
    OpenFile(&ofp, outfile, "w")
    mfp = NULL
    if (!Keyword(maskfile, "none"))
        OpenFile(&mfp, maskfile, "r")

    AllocateFloat(&phase, xsize*ysize, "phase data")
    AllocateFloat(&soln, xsize*ysize, "unwrapped data")
    AllocateByte(&bitflags, xsize*ysize, "bitflag data")

    #    READ AND PROCESS DATA
    print("Reading phase data...")
    GetPhase(in_format, ifp, infile, phase, xsize, ysize)

    # mask data
    print("Processing mask data...")
    if (mfp) {
        ReadByte(mfp, bitflags, xsize*ysize, maskfile)

    else {
        for (k=0; k<xsize*ysize; k += 1)
            bitflags[k] = 255

    for (k=0; k<xsize*ysize; k += 1)
        bitflags[k] = (!bitflags[k]) ? BORDER : 0

    FattenMask(bitflags, BORDER, 1, xsize, ysize)

    #    LOCATE AND PROCESS RESIDUES
    # compute residues and store in bitflags array
    NumRes = Residues(phase, bitflags, POS_RES, NEG_RES,
                                        BORDER, xsize, ysize)
    print("%d Residues", NumRes)

    if (debug_flag) {
        char filename[300]
        sprintf(filename, "%s.res", outfile)
        SaveByteToImage(bitflags, "residues", filename, xsize, ysize,
                                        1, 1, 0)


    #    GENERATE BRANCH CUTS
    if (dipole_flag) {    # elimate dipole-residues first
        Dipole(bitflags, xsize, ysize, BRANCH_CUT)
        i = Residues(phase, bitflags, 0, 0, BORDER | BRANCH_CUT,
                                 xsize, ysize)
        print("%d Residues are left", i)


    if (MaxCutLen==0) MaxCutLen = (xsize + ysize)/2
    GoldsteinBranchCuts(bitflags, MaxCutLen, NumRes, xsize, ysize,
                                            BRANCH_CUT)
    if (debug_flag) {
        char filename[300]
        sprintf(filename, "%s.brc", outfile)
        SaveByteToImage(bitflags, "branch cuts", filename, 
                                        xsize, ysize, 1, 1, BRANCH_CUT | BORDER)


    #    UNWRAP AROUND CUTS
    for (k=0, i=0; i<xsize*ysize; i += 1) 
        if (bitflags[i] & BRANCH_CUT) k += 1
    print("%d BRANCH CUT PIXELS", k)
    print("Unwrapping around branch cuts")
    k = UnwrapAroundCuts(phase, bitflags, soln, xsize, ysize, AVOID,
                                             0, NULL)
    if (k > 1) print("%d disconnected pieces.", k)
    else print("%d piece", k)
    print("\nFinished")
    PrintMinAndMax(xsize, ysize, soln, "solution")

    #    SAVE RESULT
    for (k=0; k<xsize*ysize; k += 1)
        soln[k] *= TWOPI; # scale output
    print("Saving unwrapped surface to file '%s'", outfile)
    WriteFloat(ofp, soln, xsize*ysize, outfile)
    delete(soln)
    delete(phase)
    delete(bitflags)

#
 * mainjump.c - Generate image of discontinuities 
 *
 * Source code files required:
 *            maskfat.c        mainjump.c                util.c

#include <stdio.h>
#include <math.h>
#include "file.h"
#include "pi.h"
#include "util.h"
#include "maskfat.h"

main (int argc, char *argv[])
{
    int                        i, j, k, m, n, a, b
    int                        ii, jj, kk
    FILE                     *ifp, *ofp, *mfp=NULL
    float                    *surf;         # array
    unsigned char    *mask=NULL, *out
    char                     infile[200], outfile[200]
    char                     maskfile[200]
    int                        xsize, ysize;     # dimensions of arrays
    double                 r, r1, r2, sum, xsum

    print("Generate discontinuity map of surface data")
    if (argc!=6) {
        fprintf(stderr,"Usage: %s surf image xsize ysize mask",
                                    argv[0])
        fprintf(stderr,"\nwhere surf = input surface (float) file")
        fprintf(stderr,"                image = output image (byte) file")
        fprintf(stderr,"                xsize, ysize = array dimensions")
        fprintf(stderr,"                mask = mask file (0=mask).")
        fprintf(stderr,"                (If no mask, enter 'none'.)")
        exit(BAD_USAGE)

    sscanf(argv[1], "%s", infile)
    sscanf(argv[2], "%s", outfile)
    sscanf(argv[3], "%d", &xsize)
    sscanf(argv[4], "%d", &ysize)
    sscanf(argv[5], "%s", maskfile)

    print("Input file = %s", infile)
    print("Output file = %s", outfile)
    print("Mask file = %s", maskfile)
    print("Array size = %d x %d", xsize, ysize)
    if (Keyword(maskfile, "none")) print("No mask file")

    #    OPEN FILES, ALLOCATE MEMORY
    OpenFile(&ifp, infile, "r"); 
    OpenFile(&ofp, outfile, "w"); 
    if (!Keyword(maskfile, "none")) OpenFile(&mfp, maskfile, "r"); 
    AllocateFloat(&surf, xsize*ysize, "surface data")
    AllocateByte(&out, xsize*ysize, "binary image")
    if (mfp) AllocateByte(&mask, xsize*ysize, "mask data")
 
    #    READ DATA
    print("Reading data...")
    ReadFloat(ifp, surf, xsize*ysize, infile)
    if (mfp) ReadByte(mfp, mask, xsize*ysize, maskfile)

    # FATTEN MASK BY 1 PIXEL (OPTIONAL)
    if (mask) {
        for (k=0; k<xsize*ysize; k += 1) mask[k] = !mask[k]
        FattenMask(mask, 1, 1, xsize, ysize)
        for (k=0; k<xsize*ysize; k += 1) mask[k] = !mask[k]


    # Generate image of discontinuities (where derivs > PI)
    for (j=0; j<ysize-1; j += 1) {
        for (i=0; i<xsize-1; i += 1) {
            k = j*xsize + i
            out[k] = 255; 
            if (mask and (!mask[k] or !mask[k+1] or !mask[k+xsize])) {
                out[k] = 0

            else {
                r = surf[k] - surf[k+1]
                if (r > PI or r < -PI) out[k] = 0
                r = surf[k] - surf[k+xsize]
                if (r > PI or r < -PI) out[k] = 0
    

 
    # Compute measure (sum) of discontinuities
    for (j=0, sum=xsum=0.0; j<ysize-1; j += 1) {
        for (i=0; i<xsize-1; i += 1) {
            k = j*xsize + i
            if (mask and (!mask[k] or !mask[k+1] or !mask[k+xsize]))
                continue
            if (!out[k])  += 1sum
             += 1xsum


    print("Saving image to %s", outfile)
    print("L0-measure (normalized to 0-1) of discontinuities %lf",
                 sum/xsum)
    print("(%lf%% of %d nonmasked pixels)", 
                 100.0*sum/xsum, (int)xsum); 
    WriteByte(ofp, out, xsize*ysize, outfile)
    delete(surf)
    fclose(ifp)
    delete(out)
    fclose(ofp)
    if (mfp) delete(mask)
    if (mfp) fclose(mfp)


#
 * mainlpno.c -- phase unwrapping by means of min. Lp norm algorithm
 *
 * Source code files required:
 *         congruen.c                 dct.c        dxdygrad.c         extract.c
 *            getqual.c                grad.c             histo.c         laplace.c
 *             lpnorm.c        mainlpno.c         maskfat.c                 pcg.c
 *         qualgrad.c        qualpseu.c         qualvar.c            raster.c
 *         residues.c         solncos.c                util.c

#include <stdio.h>
#include <math.h>
#include "file.h"
#include "histo.h"
#include "pi.h"
#include "getqual.h"
#include "congruen.h"
#include "maskfat.h"
#include "util.h"
#include "extract.h"
#include "getqual.h"
#include "lpnorm.h"
#define BORDER     0x20
#define DEFAULT_NUM_ITER    10
#define DEFAULT_PCG_ITER    20
#define DEFAULT_E0                0.001

main (int argc, char *argv[])
{
    int                        i, j, k, n
    FILE                     *ifp, *ofp, *mfp=0, *qfp=0
    float                    *phase;         # array
    float                    *soln;            # array
    float                    *qual_map;    # array
    float                    *rarray;        # array
    float                    *zarray;        # array
    float                    *parray;        # array
    float                    *dxwts;         # array
    float                    *dywts;         # array
    unsigned char    *bitflags
    double                 *xcos, *ycos
    char                     buffer[200], tempstr[200]
    char                     infile[200], outfile[200]
    char                     bmaskfile[200], qualfile[200]
    char                     format[200], modekey[200]
    int                        in_format, debug_flag
    int                        xsize, ysize;     # dimensions of arrays
    int                        xsize_actual, ysize_actual, xsize_dct, ysize_dct
    int                        avoid_code, thresh_flag, fatten, tsize
    int                        num_iter=DEFAULT_NUM_ITER
    int                        pcg_iter=DEFAULT_PCG_ITER
    double                 rmin, rmax, rscale, e0=DEFAULT_E0
    double                 one_over_twopi = 1.0/TWOPI
    UnwrapMode         mode
    char                     use[] =         # define usage statement
        "Usage: prog-name -input file -format fkey -output file"
        "    -xsize x -ysize y [ -mode mkey -bmask file -corr file"
        "    -tsize size -debug yes/no -iter num -pcg_iter nump"
        "    -e0 e0val -thresh yes/no -fat n ]"
        "where 'fkey' is a keyword designating the input file type"
        "(key = complex8, complex4, float or byte), 'x' and 'y' are"
        "the dimensions of the file, bmask is an optional byte-file"
        "for masking out undefined phase, corr is an optional byte-"
        "file of cross-correlation values, tsize is the size of the"
        "square template for averaging the corr file or quality"
        "values (default = 1), and 'mkey' is a keyword designating"
        "the type of quality map.    The value of mkey may be"
        "'min_grad' for Minimum Gradient, 'min_var' for Minimum"
        "Variance, 'max_corr' for Maximum Correlation, 'max_pseu' for"
        "Maximum Pseudocorrelation, or 'none' for no quality map." 
        "All files are simple raster files, and the output file"
        "is a raster file of elevation values of the unwrapped"
        "surface.    If the 'debug' parm is 'yes', then intermediate"
        "byte-files are saved.    The maximum number of iterations is"
        "'num', and the number of PCG iterations is 'nump'.    The"
        "normalization parameter is 'e0'.    To apply an automatic"
        "threshold to the quality map, the 'thresh' parm should be"
        "yes.    To thicken the quality mask, the 'fat' parm should be"
        "the number of pixels by which to fatten it."
    
    print("Phase Unwrapping by Minimum Lp Norm Algorithm")
    
    # GET COMMAND LINE PARAMETERS AND CHECK
    CommandLineParm(argc,argv, "-input", StringParm, infile, 1, use)
    CommandLineParm(argc,argv, "-format", StringParm, format, 1, use)
    CommandLineParm(argc,argv, "-output", StringParm, outfile, 1,use)
    CommandLineParm(argc,argv, "-xsize", IntegerParm, &xsize, 1, use)
    CommandLineParm(argc,argv, "-ysize", IntegerParm, &ysize, 1, use)
    if (!CommandLineParm(argc, argv, "-mode", StringParm,
                modekey, 0, use)) strcpy(modekey, "none")
    if (!CommandLineParm(argc, argv, "-corr", StringParm,
                qualfile, 0, use)) strcpy(qualfile, "none")
    if (!CommandLineParm(argc, argv, "-tsize", IntegerParm,
                &tsize, 0, use)) tsize = 1
    if (!CommandLineParm(argc, argv, "-bmask", StringParm, bmaskfile,
                0, use)) strcpy(bmaskfile, "none")
    if (!CommandLineParm(argc, argv, "-debug", StringParm, tempstr,
                0, use)) debug_flag = 0
    else debug_flag = Keyword(tempstr, "yes")
    if (!CommandLineParm(argc, argv, "-thresh", StringParm,
                tempstr, 0, use)) thresh_flag = 0
    else thresh_flag = Keyword(tempstr, "yes")
    if (!CommandLineParm(argc, argv, "-fat", IntegerParm,
                &fatten, 0, use)) fatten = 0
    CommandLineParm(argc, argv, "-iter", IntegerParm, &num_iter,
                                    0, use)
    if (num_iter < 0) num_iter = 1
    CommandLineParm(argc, argv, "-pcg_iter", IntegerParm, &pcg_iter,
                                    0, use)
    if (pcg_iter < 0) pcg_iter = 1
    CommandLineParm(argc, argv, "-e0", DoubleParm, &e0, 0, use)

    if (Keyword(format, "complex8"))    in_format = 0
    else if (Keyword(format, "complex4"))    in_format = 1
    else if (Keyword(format, "byte"))    in_format = 2
    else if (Keyword(format, "float"))    in_format = 3
    else {
        fprintf(stderr, "Unrecognized format: %s", format)
        exit(BAD_PARAMETER)


    print("Input file =    %s", infile)
    print("Input file type = %s", format)
    print("Output file =    %s", outfile)
    print("File dimensions = %dx%d (cols x rows).", xsize, ysize)
    if (Keyword(bmaskfile, "none")) print("No border mask file.")
    else print("Border mask file = %s", bmaskfile)
    print("Quality mode = %s", modekey)
    if (modekey==none) {
        print("No quality map.")
        strcpy(qualfile, "none")

    else {
        if (Keyword(qualfile, "none")) print("No correlation file.")
        else print("Correlation image file = %s", qualfile)
        print("Averaging template size = %d", tsize)
        if (tsize < 0 or tsize > 30) {
            fprintf(stderr, "Illegal size: must be between 0 and 30")
            exit(BAD_PARAMETER)


    mode = SetQualityMode(modekey, qualfile, 1)
    if (mode < 0) exit(BAD_PARAMETER);    # error msg already printed

    # Increase dimensions to power of two (plus one)
    xsize_actual = xsize
    ysize_actual = ysize
    for (xsize_dct = 1; xsize_dct+1 < xsize_actual; xsize_dct*=2)
        
    xsize_dct+=1
    for (ysize_dct = 1; ysize_dct+1 < ysize_actual; ysize_dct*=2)
        
    ysize_dct+=1
    if (xsize_dct != xsize_actual or ysize_dct != ysize_actual) {
        print("Dimensions increase from %dx%d to %dx%d for FFT's",
                     xsize_actual, ysize_actual, xsize_dct, ysize_dct)


    #    OPEN FILES, ALLOCATE MEMORY
    # note: xsize_dct >= xsize;    ysize_dct >= ysize
    xsize = xsize_dct
    ysize = ysize_dct
    AllocateFloat(&phase, xsize*ysize, "phase data")
    AllocateFloat(&soln, xsize*ysize, "scratch data")
    AllocateFloat(&qual_map, xsize*ysize, "quality map")
    AllocateByte(&bitflags, xsize*ysize, "bitflags array")
    OpenFile(&ifp, infile, "r")
    OpenFile(&ofp, outfile, "w")
    if (!Keyword(bmaskfile, "none")) OpenFile(&mfp, bmaskfile, "r")
    if (mode==corr_coeffs) OpenFile(&qfp, qualfile, "r")

    #    READ AND PROCESS DATA
    xsize = xsize_actual
    ysize = ysize_actual
    print("Reading input data...")
    GetPhase(in_format, ifp, infile, phase, xsize, ysize)

    if (qfp) {
        print("Reading quality data...")
        # borrow the bitflags array temporarily
        ReadByte(qfp, bitflags, xsize*ysize, qualfile)
        # process data and store in quality map array
        AverageByteToFloat(bitflags, qual_map, tsize, xsize, ysize)


    # border mask data
    print("Processing border mask data...")
    if (mfp) {
        ReadByte(mfp, bitflags, xsize*ysize, bmaskfile)

    else {
        for (k=0; k<xsize*ysize; k += 1)
            bitflags[k] = 255

    for (k=0; k<xsize*ysize; k += 1) {
        bitflags[k] = (!bitflags[k]) ? BORDER : 0

    if (mfp) FattenMask(bitflags, BORDER, 1, xsize, ysize)

    GetQualityMap(mode, qual_map, phase, bitflags, BORDER,
                                tsize, xsize, ysize)
    if (thresh_flag) {
        HistoAndThresh(qual_map, xsize, ysize, 0.0, 0, 0.0,
                                     bitflags, BORDER)
        if (fatten > 0)
            FattenQual(qual_map, fatten, xsize, ysize)

    # Set border (masked) pixels to 0-wts
    for (k=0; k<xsize*ysize; k += 1) {
        if (bitflags[k] & BORDER) qual_map[k] = 0.0


    if (debug_flag) {
        char filename[300]
        sprintf(filename, "%s.qual", outfile)
        SaveFloatToImage(qual_map, "quality", filename, xsize, ysize,
                                         0, 0, 0)


    # embed arrays in possibly larger FFT/DCT arrays
    for (j=ysize_dct-1; j>=0; j--) {
        for (i=xsize_dct-1; i>=0; i--) {
            if (i<xsize_actual and j<ysize_actual)
                phase[j*xsize_dct + i] = phase[j*xsize_actual + i]
            else phase[j*xsize_dct + i] = 0.0


    for (j=ysize_dct-1; j>=0; j--) {
        for (i=xsize_dct-1; i>=0; i--) {
            if (i<xsize_actual and j<ysize_actual)
                bitflags[j*xsize_dct + i] = bitflags[j*xsize_actual + i]
            else bitflags[j*xsize_dct + i] = BORDER


    if (qual_map) {
        for (j=ysize_dct-1; j>=0; j--) {
            for (i=xsize_dct-1; i>=0; i--) {
                if (i<xsize_actual and j<ysize_actual)
                    qual_map[j*xsize_dct + i] = qual_map[j*xsize_actual + i]
                else
                    qual_map[j*xsize_dct + i] = 0.0




    # Set dimensions to DCT dimensions
    xsize = xsize_dct
    ysize = ysize_dct
    # Allocate more memory
    AllocateFloat(&rarray, xsize*ysize, "r array data")
    AllocateFloat(&zarray, xsize*ysize, "z array data")
    AllocateFloat(&parray, xsize*ysize, "p array data")
    AllocateFloat(&dxwts, xsize*ysize, "dx weights")
    AllocateFloat(&dywts, xsize*ysize, "dy weights")

    #    UNWRAP
    for (k=0; k<xsize*ysize; k += 1) soln[k] = 0.0
    print("Unwrapping...")
    LpNormUnwrap(soln, phase, dxwts, dywts, bitflags, qual_map,
        rarray, zarray, parray, num_iter, pcg_iter, e0, xsize, ysize)
    print("\nFinished")
    PrintMinAndMax(xsize, ysize, soln, "solution")

    # restore dimensions to input sizes and save solution
    xsize = xsize_actual
    ysize = ysize_actual
    for (j=0; j<ysize; j += 1) {
        for (i=0; i<xsize; i += 1) {
            soln[j*xsize + i] = TWOPI*soln[j*xsize_dct + i]


    # save solution
    print("Saving unwrapped surface to file '%s'", outfile)
    WriteFloat(ofp, soln, xsize*ysize, outfile)
    delete(rarray)
    delete(zarray)
    delete(parray)
    delete(soln)
    delete(qual_map)
    delete(bitflags)
    delete(dxwts)
    delete(dywts)

#
 * mainmcut.c - phase unwrapping with quality-guided mask cuts
 *
 * Source code files required:
 *         dxdygrad.c         extract.c         getqual.c                grad.c
 *                 list.c        mainmcut.c         maskcut.c         maskfat.c
 *         maskthin.c                path.c        qualgrad.c        qualpseu.c
 *            qualvar.c        residues.c                util.c

#include <stdio.h>
#include <math.h>
#include "file.h"
#include "pi.h"
#include "list.h"
#include "util.h"
#include "extract.h"
#include "getqual.h"
#include "maskcut.h"
#include "residues.h"
#include "path.h"

main (int argc, char *argv[])
{
    int                        i, j, k, NumRes, mask_code
    FILE                     *ifp, *ofp, *mfp=0, *qfp=0
    float                    *phase;         # array
    float                    *soln;            # array
    float                    *qual_map;    # array
    unsigned char    *bitflags;    # array
    char                     buffer[200], string[200], infile[200]
    char                     outfile[200], bmaskfile[200], qualfile[200]
    char                     tempstr[200], format[200], modekey[200]
    int                        xsize, ysize;     # dimensions of arrays
    int                        in_format, debug_flag, avoid_code, tsize
    double                 minqual, maxqual, rmin, rmax, rscale
    UnwrapMode         mode
    char                     use[]    =        # define usage statement
        "Usage: prog-name -input file -format fkey -output file"
        "    -xsize x -ysize y -mode mkey [ -bmask file -corr file"
        "    -tsize size -debug yes/no]"
        "where 'fkey' is a keyword designating the input file type"
        "(key = complex8, complex4, float or byte), 'x' and 'y' are"
        "the dimensions of the file, bmask is an optional byte-file"
        "of masks for masking out undefined phase values, corr"
        "is an optional byte-file of cross-correlation values,"
        "tsize is the size of the square template for averaging the"
        "corr file or quality values (default = 1), and 'mkey' is"
        "a keyword designating the source of the quality values"
        "that guide the masking path.    The value of mkey must be"
        "'min_grad' for minimum gradients, 'min_var' for minimum"
        "phase variances, 'max_corr' for maximum cross"
        "correlations (i.e., corr file), or 'max_pseu' for max-"
        "imum pseudocorrelations.    All of the files are simple"
        "raster files, and the output file consists of floating"
        "point numbers that define the heights of the unwrapped"
        "surface.    If the 'debug' parm is 'yes', then intermediate"
        "byte-files are saved (residues, unwrapping paths, etc.)"
    
    print("Phase Unwrapping By Quality-Guided Mask Cuts")

    # GET COMMAND LINE PARAMETERS AND CHECK
    CommandLineParm(argc,argv, "-input", StringParm, infile, 1, use); 
    CommandLineParm(argc,argv, "-format", StringParm, format, 1, use); 
    CommandLineParm(argc,argv, "-output", StringParm, outfile, 1,use); 
    CommandLineParm(argc,argv, "-xsize", IntegerParm, &xsize, 1, use); 
    CommandLineParm(argc,argv, "-ysize", IntegerParm, &ysize, 1, use); 
    CommandLineParm(argc,argv, "-mode", StringParm, modekey, 1, use); 
    if (!CommandLineParm(argc, argv, "-bmask", StringParm, bmaskfile,
                0, use)) strcpy(bmaskfile, "none")
    if (!CommandLineParm(argc, argv, "-corr", StringParm, qualfile,
                0, use)) strcpy(qualfile, "none")
    if (!CommandLineParm(argc, argv, "-tsize", IntegerParm, &tsize,
                0, use)) tsize = 1
    if (!CommandLineParm(argc, argv, "-debug", StringParm, tempstr,
                0, use)) debug_flag = 0
    else debug_flag = Keyword(tempstr, "yes")

    if (Keyword(format, "complex8"))    in_format = 0
    else if (Keyword(format, "complex4"))    in_format = 1
    else if (Keyword(format, "byte"))    in_format = 2
    else if (Keyword(format, "float"))    in_format = 3
    else {
        fprintf(stderr, "Unrecognized format: %s", format)
        exit(BAD_PARAMETER)


    print("Input file =    %s", infile)
    print("Input file type = %s", format); 
    print("Output file =    %s", outfile)
    print("File dimensions = %dx%d (cols x rows).", xsize, ysize)

    if (Keyword(bmaskfile, "none")) print("No border mask file.")
    else print("Border mask file = %s", bmaskfile)
    if (Keyword(qualfile, "none")) print("No quality file.")
    else print("Correlation image file = %s", qualfile)

    print("Averaging template size = %d", tsize)
    if (tsize < 0 or tsize > 10) { 
        fprintf(stderr, "Illegal size: must be between 0 and 10")
        exit(BAD_PARAMETER)

    print("Quality mode = %s", modekey)
    mode = SetQualityMode(modekey, qualfile, 1)
    if (mode < 0) exit(BAD_PARAMETER);    # error msg already printed

    #    OPEN FILES, ALLOCATE MEMORY
    OpenFile(&ifp, infile, "r")
    OpenFile(&ofp, outfile, "w")
    if (!Keyword(bmaskfile, "none"))
        OpenFile(&mfp, bmaskfile, "r")
    if (mode==corr_coeffs) 
        OpenFile(&qfp, qualfile, "r")

    AllocateFloat(&phase, xsize*ysize, "phase data")
    AllocateFloat(&soln, xsize*ysize, "unwrapped data")
    AllocateByte(&bitflags, xsize*ysize, "bitflags array")
    AllocateFloat(&qual_map, xsize*ysize, "quality map")

    #    READ AND PROCESS DATA
    print("Reading phase data...")
    GetPhase(in_format, ifp, infile, phase, xsize, ysize)

    if (qfp) {
        print("Reading quality data...")
        # borrow the bitflags array temporarily
        ReadByte(qfp, bitflags, xsize*ysize, qualfile)
        # process data and store in qual_map
        AverageByteToFloat(bitflags, qual_map, tsize, xsize, ysize)


    # border mask data
    print("Processing border mask data...")
    if (mfp) {
        ReadByte(mfp, bitflags, xsize*ysize, bmaskfile)

    else { 
        for (k=0; k<xsize*ysize; k += 1)
            bitflags[k] = 255

    for (k=0; k<xsize*ysize; k += 1)
        bitflags[k] = (!bitflags[k]) ? BORDER : 0
    FattenMask(bitflags, BORDER, 1, xsize, ysize)

    GetQualityMap(mode, qual_map, phase, bitflags, BORDER,
                                tsize, xsize, ysize)

    # initialize soln array
    for (k=0; k<xsize*ysize; k += 1) 
        soln[k] = 0

    #    LOCATE AND PROCESS RESIDUES
    # compute residues and store in bitflags array
    NumRes = Residues(phase, bitflags, POS_RES, NEG_RES, BORDER,
                                        xsize, ysize)
    print("%d Residues", NumRes)

    if (debug_flag) {
        char filename[300]
        sprintf(filename, "%s.res", outfile); 
        SaveByteToImage(bitflags, "residues", filename, xsize, ysize,
                                        1, 1, 0)


    # mark residues and border in quality map array
    for (k=0, maxqual = - (minqual = 1.0E+10); k<xsize*ysize; k += 1) {
        if (maxqual < qual_map[k]) maxqual = qual_map[k]
        if (minqual > qual_map[k]) minqual = qual_map[k]

    minqual -= 0.0001
    maxqual += 0.0001
    # change quality values to cost values
    # (to follow low-quality mask cuts)
    # i.e., reverse high and low qualities
    for (k=0; k<xsize*ysize; k += 1) 
        qual_map[k] = maxqual - qual_map[k]
    # mark residues and border as high quality
    for (j=0; j<ysize; j += 1) {
        for (i=0; i<xsize; i += 1) {
            k = j*xsize + i
            if (i==0 or i==xsize-1 or j==0 or j==ysize-1) 
                qual_map[k] = maxqual
            else if (bitflags[k] & BORDER) 
                qual_map[k] = maxqual
            else if (bitflags[k] & RESIDUE) 
                qual_map[k] = maxqual



    if (debug_flag) {
        char filename[300]
        sprintf(filename, "%s.qual", outfile); 
        SaveFloatToImage(qual_map, "quality map", filename,
                                         xsize, ysize, 1, 0, 0)


    #    GENERATE MASK CUTS
    if (mode==gradient and tsize==1) mode = dxdygrad
    mask_code = BRANCH_CUT
    k = QualityGuidedMask(phase, qual_map, bitflags, xsize, ysize,
                                                mask_code, mode, debug_flag, outfile)
    if (k > 1) print("\n%d pieces", k)
    if (debug_flag) {
        char filename[300]
        sprintf(filename, "%s.fatmask", outfile); 
        SaveByteToImage(bitflags, "fat masks", filename, xsize, ysize,
                                        1, 1, mask_code | BORDER)

    for (k=0; k<xsize*ysize; k += 1) {
        bitflags[k] &= (mask_code | BORDER | RESIDUE)


    #    THIN MASK CUTS
    ThinMask(bitflags, xsize, ysize, mask_code)

    if (debug_flag) {
        char filename[300]
        sprintf(filename, "%s.thinmask", outfile); 
        SaveByteToImage(bitflags, "thin masks", filename, xsize, ysize,
                                        1, 1, mask_code | BORDER)


    #    UNWRAP AROUND CUTS
    for (k=0; k<xsize*ysize; k += 1)    # initialize solution array
        soln[k] = 0.0
    print("Unwrapping around branch cuts")
    k = UnwrapAroundCuts(phase, bitflags, soln, xsize, ysize,
                                             mask_code, debug_flag, outfile)
    if (k > 1) print("%d disconnected pieces.", k)
    else print("%d piece", k)
    
    print("\nFinished")
    PrintMinAndMax(xsize, ysize, soln, "solution")

    #    SAVE RESULT
    # scale output
    for (k=0; k<xsize*ysize; k += 1) 
        soln[k] *= TWOPI
    print("Saving unwrapped surface to file '%s'", outfile)
    WriteFloat(ofp, soln, xsize*ysize, outfile)
    delete(soln)
    delete(phase)
    delete(bitflags)
    delete(qual_map)

#
 * mainpcg.c -- phase unwrapping by means of PCG algorithm
 *
 * Source code files required:
 *         congruen.c                 dct.c        dxdygrad.c         extract.c
 *            getqual.c                grad.c             histo.c         laplace.c
 *            mainpcg.c         maskfat.c                 pcg.c        qualgrad.c
 *         qualpseu.c         qualvar.c         solncos.c                util.c

#include <stdio.h>
#include <math.h>
#include "file.h"
#include "histo.h"
#include "pi.h"
#include "getqual.h"
#include "congruen.h"
#include "maskfat.h"
#include "util.h"
#include "extract.h"
#include "pcg.h"
#include "laplace.h"
#define BORDER     0x20
#define DEFAULT_NUM_ITER    20

main (int argc, char *argv[])
{
    int                        i, j, k, n
    FILE                     *ifp, *ofp, *mfp=0, *qfp=0
    float                    *phase;         # array
    float                    *soln;            # array
    float                    *qual_map;    # array
    float                    *rarray;    # array
    float                    *parray;    # array
    float                    *zarray;    # array
    unsigned char    *bitflags
    double                 *xcos, *ycos
    char                     buffer[200], tempstr[200]
    char                     infile[200], outfile[200]
    char                     bmaskfile[200], qualfile[200]
    char                     format[200], modekey[200]
    int                        in_format, debug_flag
    int                        xsize, ysize;     # dimensions of arrays
    int                        xsize_actual, ysize_actual
    int                        xsize_dct, ysize_dct
    int                        avoid_code, thresh_flag, fatten
    int                        tsize, num_iter=DEFAULT_NUM_ITER
    double                 rmin, rmax, rscale, epsi_con
    double                 one_over_twopi = 1.0/TWOPI
    UnwrapMode         mode
    char                     use[] =             # define usage statement
     "Usage: prog-name -input file -format fkey -output file"
     "    -xsize x -ysize y [ -mode mkey -bmask file -corr file"
     "    -tsize size -debug yes/no -iter num -converge epsi"
     "    -thresh yes/no -fat n ]"
     "where 'fkey' is a keyword designating the input file type"
     "(key = complex8, complex4, float or byte), 'x' and 'y' are"
     "the dimensions of the file, bmask is an optional byte-file"
     "of masks for masking out undefined phase values, corr"
     "is an optional byte-file of cross-correlation values,"
     "tsize is the size of the square template for averaging the"
     "corr file or quality values (default = 1), and 'mkey' is"
     "a keyword designating the source of the quality values"
     "that guide the unwrapping path.    The value of mkey may be"
     "'min_grad' for Minimum Gradient unwrapping, 'min_var' for"
     "Minimum Variance unwrapping, 'max_corr' for Maximum"
     "Correlation unwrapping, or 'max_pseu' for Maximum"
     "Pseudocorrelation unwrapping.    All files are simple"
     "raster files, and the output file consists of floating"
     "point numbers that define the heights of the unwrapped"
     "surface.    If the 'debug' parm is 'yes', then intermediate"
     "byte-files are saved (quality map, etc.)    The maximum number"
     "of iterations is given by 'num' (default 20), and the con-"
     "vergence tolerance is given by 'epsi' (default = 0).    To"
     "apply an automatic threshold to the quality map to make"
     "a quality mask, the 'thresh' parm should be yes.    To"
     "thicken the quality mask, the 'fat' parm should be the"
     "number of pixels by which to thicken."

    print("Phase Unwrapping by PCG algorithm")

    # GET COMMAND LINE PARAMETERS AND CHECK
    CommandLineParm(argc,argv, "-input", StringParm, infile, 1, use)
    CommandLineParm(argc,argv, "-format", StringParm, format, 1, use)
    CommandLineParm(argc,argv, "-output", StringParm, outfile, 1,use)
    CommandLineParm(argc,argv, "-xsize", IntegerParm, &xsize, 1, use)
    CommandLineParm(argc,argv, "-ysize", IntegerParm, &ysize, 1, use)
    if (!CommandLineParm(argc, argv, "-mode", StringParm,
                modekey, 0, use)) strcpy(modekey, "none")
    if (!CommandLineParm(argc, argv, "-bmask", StringParm,        
                bmaskfile, 0, use)) strcpy(bmaskfile, "none")
    if (!CommandLineParm(argc, argv, "-corr", StringParm,
                qualfile, 0, use)) strcpy(qualfile, "none")
    if (!CommandLineParm(argc, argv, "-tsize", IntegerParm,
                &tsize, 0, use)) tsize = 1
    if (!CommandLineParm(argc, argv, "-debug", StringParm,
                tempstr, 0, use)) debug_flag = 0
    else debug_flag = Keyword(tempstr, "yes")
    CommandLineParm(argc, argv, "-iter", IntegerParm,
                                    &num_iter, 0, use)
    if (num_iter < 0) num_iter = 1
    if (!CommandLineParm(argc, argv, "-toler", DoubleParm,
                &epsi_con, 0, use)) epsi_con = 0.0
    if (!CommandLineParm(argc, argv, "-thresh", StringParm,
                tempstr, 0, use)) thresh_flag = 0
    else thresh_flag = Keyword(tempstr, "yes")
    if (!CommandLineParm(argc, argv, "-fat", IntegerParm,
                &fatten, 0, use)) fatten = 0

    if (Keyword(format, "complex8"))    in_format = 0
    else if (Keyword(format, "complex4"))    in_format = 1
    else if (Keyword(format, "byte"))    in_format = 2
    else if (Keyword(format, "float"))    in_format = 3
    else {
        fprintf(stderr, "Unrecognized format: %s", format)
        exit(BAD_PARAMETER)


    print("Input file =    %s", infile)
    print("Input file type = %s", format)
    print("Output file =    %s", outfile)
    print("File dimensions = %dx%d (cols x rows).", xsize, ysize)

    if (Keyword(bmaskfile, "none")) print("No border mask file.")
    else print("Border mask file = %s", bmaskfile)
    if (Keyword(qualfile, "none")) print("No quality file.")
    else print("Correlation image file = %s", qualfile)

    print("Averaging template size = %d", tsize)
    if (tsize < 0 or tsize > 30) {
        fprintf(stderr, "Illegal size: must be between 0 and 30")
        exit(BAD_PARAMETER)


    print("Quality mode = %s", modekey)
    mode = SetQualityMode(modekey, qualfile, 0)
    if (mode < 0) exit(BAD_PARAMETER);    # error msg already printed

    # Increase dimensions to power of two (plus one)
    xsize_actual = xsize
    ysize_actual = ysize
    for (xsize_dct = 1; xsize_dct + 1 < xsize_actual; xsize_dct*=2)
        
    xsize_dct = xsize_dct + 1
    for (ysize_dct = 1; ysize_dct + 1 < ysize_actual; ysize_dct*=2)
        
    ysize_dct = ysize_dct + 1
    if (xsize_dct != xsize_actual or ysize_dct != ysize_actual) {
        print("Dim's increased from %dx%d to %dx%d for FFT/DCT's",
                     xsize_actual, ysize_actual, xsize_dct, ysize_dct)


    #    OPEN FILES, ALLOCATE MEMORY
    # note: xsize_dct >= xsize;    ysize_dct >= ysize
    xsize = xsize_dct
    ysize = ysize_dct
    AllocateFloat(&phase, xsize*ysize, "phase data")
    AllocateFloat(&soln, xsize*ysize, "scratch data")
    AllocateFloat(&qual_map, xsize*ysize, "quality map")
    AllocateByte(&bitflags, xsize*ysize, "bitflags array")
    OpenFile(&ifp, infile, "r")
    OpenFile(&ofp, outfile, "w")
    if (!Keyword(bmaskfile, "none")) OpenFile(&mfp, bmaskfile, "r")
    if (mode==corr_coeffs) OpenFile(&qfp, qualfile, "r")

    #    READ AND PROCESS DATA
    xsize = xsize_actual
    ysize = ysize_actual
    print("Reading input data...")
    GetPhase(in_format, ifp, infile, phase, xsize, ysize)

    if (qfp) {
        print("Reading quality data...")
        # borrow the bitflags array temporarily
        ReadByte(qfp, bitflags, xsize*ysize, qualfile)
        # process data and store in quality map array
        AverageByteToFloat(bitflags, qual_map, tsize, xsize, ysize)


    # border mask data
    print("Processing border mask data...")
    if (mfp) {
        ReadByte(mfp, bitflags, xsize*ysize, bmaskfile)

    else {
        for (k=0; k<xsize*ysize; k += 1)
            bitflags[k] = 255

    for (k=0; k<xsize*ysize; k += 1) {
        bitflags[k] = (!bitflags[k]) ? BORDER : 0

    if (mfp) FattenMask(bitflags, BORDER, 1, xsize, ysize)

    GetQualityMap(mode, qual_map, phase, bitflags, BORDER,
                                tsize, xsize, ysize)
    
    if (thresh_flag) {
        HistoAndThresh(qual_map, xsize, ysize, 0.0, 0, 0.0,
                                     bitflags, BORDER)
        if (fatten > 0)
            FattenQual(qual_map, fatten, xsize, ysize)

    # Set border (masked) pixels to 0-wts and free bitflags
    for (k=0; k<xsize*ysize; k += 1) {
        if (bitflags[k] & BORDER) qual_map[k] = 0.0

    delete(bitflags);     # bitflags array no longer needed

    if (debug_flag) {
        char filename[300]
        sprintf(filename, "%s.qual", outfile)
        SaveFloatToImage(qual_map, "quality", filename, xsize, ysize,
                                         0, 0, 0)


    # embed arrays in possibly larger FFT/DCT arrays
    for (j=ysize_dct-1; j>=0; j--) {
        for (i=xsize_dct-1; i>=0; i--) {
            if (i<xsize_actual and j<ysize_actual)
                phase[j*xsize_dct + i] = phase[j*xsize_actual + i]
            else phase[j*xsize_dct + i] = 0.0


    if (qual_map) {
        for (j=ysize_dct-1; j>=0; j--) {
            for (i=xsize_dct-1; i>=0; i--) {
                if (i<xsize_actual and j<ysize_actual)
                    qual_map[j*xsize_dct + i] = qual_map[j*xsize_actual + i]
                else
                    qual_map[j*xsize_dct + i] = 0.0




    # Set dimensions to DCT dimensions
    xsize = xsize_dct
    ysize = ysize_dct
    # Allocate more memory
    AllocateFloat(&rarray, xsize*ysize, "r array data")
    AllocateFloat(&zarray, xsize*ysize, "z array data")
    AllocateFloat(&parray, xsize*ysize, "p array data")

    #    UNWRAP
    print("Unwrapping...")
    for (k=0; k<xsize*ysize; k += 1)    soln[k] = 0.0
    ComputeLaplacian(phase, rarray, qual_map, NULL, xsize, ysize, 1)
    PCGUnwrap(rarray, zarray, parray, soln, qual_map, NULL,
                        xsize, ysize, num_iter, epsi_con)
    print("\nFinished")
    PrintMinAndMax(xsize, ysize, soln, "solution")

    # restore dimensions to input sizes
    xsize = xsize_actual
    ysize = ysize_actual
    # extract results from enlarged arrays
    for (j=0; j<ysize; j += 1) {
        for (i=0; i<xsize; i += 1) {
            soln[j*xsize + i] = soln[j*xsize_dct + i]


    for (j=0; j<ysize; j += 1) {
        for (i=0; i<xsize; i += 1) {
            qual_map[j*xsize + i] = qual_map[j*xsize_dct + i]


    for (j=0; j<ysize; j += 1) {
        for (i=0; i<xsize; i += 1) {
            phase[j*xsize + i] = phase[j*xsize_dct + i]


    # make result congruent to input phase
#**    CongruentSoln(soln, phase, qual_map, xsize, ysize); **
    # scale and save result
    for (k=0; k<xsize*ysize; k += 1) soln[k] *= TWOPI
    WriteFloat(ofp, soln, xsize*ysize, outfile)
    print("Wrote congruent surface to file %s", outfile)
    delete(qual_map)
    delete(phase)
    delete(soln)
    delete(rarray)
    delete(parray)
    delete(zarray)

#
 * mainqual.c -- phase unwrapping by quality-guided path following
 *
 * Source code files required:
 *         dxdygrad.c         extract.c         getqual.c                grad.c
 *                 list.c        mainqual.c         maskfat.c         quality.c
 *         qualgrad.c        qualpseu.c         qualvar.c                util.c

#include <stdio.h>
#include <math.h>
#include "file.h"
#include "pi.h"
#include "quality.h"
#include "util.h"
#include "extract.h"
#include "getqual.h"

#define BORDER (0x20)
main (int argc, char *argv[])
{
    int                        i, j, k
    FILE                     *ifp, *ofp, *mfp=0, *qfp=0
    float                    *phase;         # array
    float                    *soln;            # array
    float                    *qual_map;    # array
    unsigned char    *bitflags;    # array
    char                     buffer[200], string[200]
    char                     infile[200], outfile[200]
    char                     bmaskfile[200], qualfile[200]
    char                     tempstr[200], format[200], modekey[200]
    int                        xsize, ysize;     # dimensions of arrays
    int                        tsize, in_format, debug_flag, avoid_code
    double                 rmin, rmax, rscale
    UnwrapMode         mode
    char                     use[] =            # define usage statement
        "Usage: prog-name -input file -format fkey -output file"
        "    -xsize x -ysize y -mode mkey [ -bmask file -corr file"
        "    -tsize size -debug yes/no]"
        "where 'fkey' is a keyword designating the input file type"
        "(key = complex8, complex4, float or byte), 'x' and 'y' are"
        "the dimensions of the file, bmask is an optional byte-file"
        "of masks for masking out undefined phase values, corr"
        "is an optional byte-file of cross-correlation values,"
        "tsize is the size of the square template for averaging the"
        "corr file or quality values (default = 1), and 'mkey' is"
        "a keyword designating the source of the quality values"
        "that guide the unwrapping path.    The value of mkey must be"
        "'min_grad' for Minimum Gradient unwrapping, 'min_var' for"
        "Minimum Variance unwrapping, 'max_corr' for Maximum"
        "Correlation unwrapping, or 'max_pseu' for Maximum"
        "Pseudocorrelation unwrapping.    All files are simple"
        "raster files, and the output file consists of floating"
        "point numbers that define the heights of the unwrapped"
        "surface.    If the 'debug' parm is 'yes', then intermediate"
        "byte-files are saved (quality map, unwrapping paths, etc.)"

    print("Phase Unwrapping by Quality Map Guided Path Following")

    # GET COMMAND LINE PARAMETERS AND CHECK
    CommandLineParm(argc,argv, "-input", StringParm, infile, 1, use)
    CommandLineParm(argc,argv, "-format", StringParm, format, 1, use)
    CommandLineParm(argc,argv, "-output", StringParm, outfile, 1,use)
    CommandLineParm(argc,argv, "-xsize", IntegerParm, &xsize, 1, use)
    CommandLineParm(argc,argv, "-ysize", IntegerParm, &ysize, 1, use)
    CommandLineParm(argc,argv, "-mode", StringParm, modekey, 1, use)
    if (!CommandLineParm(argc, argv, "-bmask", StringParm,
                bmaskfile, 0, use)) strcpy(bmaskfile, "none")
    if (!CommandLineParm(argc, argv, "-corr", StringParm,
                qualfile, 0, use)) strcpy(qualfile, "none")
    if (!CommandLineParm(argc, argv, "-tsize", IntegerParm,
                &tsize, 0, use)) tsize = 1
    if (!CommandLineParm(argc, argv, "-debug", StringParm,
                tempstr, 0, use)) debug_flag = 0
    else debug_flag = Keyword(tempstr, "yes")

    if (Keyword(format, "complex8"))    in_format = 0
    else if (Keyword(format, "complex4"))    in_format = 1
    else if (Keyword(format, "byte"))    in_format = 2
    else if (Keyword(format, "float"))    in_format = 3
    else {
        fprintf(stderr, "Unrecognized format: %s", format)
        exit(BAD_PARAMETER)


    print("Input file =    %s", infile)
    print("Input file type = %s", format)
    print("Output file =    %s", outfile)
    print("File dimensions = %dx%d (cols x rows).", xsize, ysize)

    if (Keyword(bmaskfile, "none")) print("No border mask file.")
    else print("Border mask file = %s", bmaskfile)
    if (Keyword(qualfile, "none")) print("No quality file.")
    else print("Correlation image file = %s", qualfile)

    print("Averaging template size = %d", tsize)
    if (tsize < 0 or tsize > 30) {
        fprintf(stderr, "Illegal size: must be between 0 and 30")
        exit(BAD_PARAMETER)


    print("Quality mode = %s", modekey)
    mode = SetQualityMode(modekey, qualfile, 0)
    if (mode < 0) exit(BAD_PARAMETER);    # error msg already printed

    #    OPEN FILES, ALLOCATE MEMORY
    OpenFile(&ifp, infile, "r")
    OpenFile(&ofp, outfile, "w")
    if (!Keyword(bmaskfile, "none"))
        OpenFile(&mfp, bmaskfile, "r")
    if (mode==corr_coeffs) 
        OpenFile(&qfp, qualfile, "r")

    AllocateFloat(&phase, xsize*ysize, "phase data")
    AllocateFloat(&soln, xsize*ysize, "unwrapped data")
    AllocateByte(&bitflags, xsize*ysize, "bitflag data")
    AllocateFloat(&qual_map, xsize*ysize, "quality map")

    #    READ AND PROCESS DATA
    print("Reading phase data...")
    GetPhase(in_format, ifp, infile, phase, xsize, ysize)

    if (qfp) {
        print("Reading quality data...")
        # borrow the bitflags array temporarily
        ReadByte(qfp, bitflags, xsize*ysize, qualfile)
        # process data and store in quality map array
        AverageByteToFloat(bitflags, qual_map, tsize, xsize, ysize)


    # border mask data
    print("Processing border mask data...")
    if (mfp) {
        ReadByte(mfp, bitflags, xsize*ysize, bmaskfile)

    else {
        for (k=0; k<xsize*ysize; k += 1)
            bitflags[k] = 255

    for (k=0; k<xsize*ysize; k += 1)
        bitflags[k] = (!bitflags[k]) ? BORDER : 0
    FattenMask(bitflags, BORDER, 1, xsize, ysize)

    GetQualityMap(mode, qual_map, phase, bitflags, BORDER,
                                tsize, xsize, ysize)
    if (debug_flag) {
        char filename[300]
        sprintf(filename, "%s.qual", outfile)
        SaveFloatToImage(qual_map, "quality", filename, xsize, ysize,
                                         0, 0, 0)


    # initialize soln array
    for (k=0; k<xsize*ysize; k += 1)
        soln[k] = 0

    #    UNWRAP
    avoid_code = BORDER
    if (mode==gradient and tsize==1) mode = dxdygrad
    k = QualityGuidedPathFollower(phase, qual_map, bitflags, soln,
                            xsize, ysize, avoid_code, mode, debug_flag, outfile)
    if (k > 1) print("\n%d pieces", k)

    print("\nFinished")
    PrintMinAndMax(xsize, ysize, soln, "solution")

    #    SAVE RESULT
    # scale output
    for (k=0; k<xsize*ysize; k += 1) 
        soln[k] *= TWOPI

    print("Saving unwrapped surface to file '%s'", outfile)
    WriteFloat(ofp, soln, xsize*ysize, outfile)
    delete(soln)
    delete(phase)
    delete(bitflags)
    delete(qual_map)

#
 * mainunmg.c -- unweighted least-squares phase unwrapping
 *                             by means of unweighted multigrid algorithm
 *
 * Source code files required:
 *         dxdygrad.c         extract.c                grad.c         gridmem.c
 *            gridops.c        mainunmg.c             relax.c             unfmg.c
 *             ungrid.c                util.c

#include <stdio.h>
#include <math.h>
#include "file.h"
#include "pi.h"
#include "util.h"
#include "extract.h"
#include "unfmg.h"
#define DEFAULT_NUM_ITER        2
#define DEFAULT_NUM_CYCLES    2

main (int argc, char *argv[])
{
    int                        i, j, k, n
    FILE                     *ifp, *ofp
    float                    *phase;         # array
    float                    *soln;            # array
    float                    *dx;                # array
    float                    *dy;                # array
    char                     infile[200], outfile[200]
    char                     format[200]
    int                        in_format
    int                        xsize, ysize;     # dimensions of arrays
    int                        num_iter=DEFAULT_NUM_ITER
    int                        num_cycles=DEFAULT_NUM_CYCLES
    char                     use[] =     # define usage statement
        "Usage: prog-name -input file -format fkey -output file"
        "    -xsize x -ysize y [ -cycles numc -iter num]"
        "where 'fkey' is a keyword designating the input file type"
        "(key = complex8, complex4, float or byte) and 'x' and 'y'"
        "are the dimensions of the file.    All files are simple"
        "raster files, and the output file consists of floating"
        "point numbers that define the heights of the unwrapped"
        "surface.    The number of cycles is given by 'numc'.    The"
        "number of Gauss-Seidel iterations is given by 'num'."

    print("Unweighted Phase Unwrapping by Multigrid Algorithm")
            
    # GET COMMAND LINE PARAMETERS AND CHECK
    CommandLineParm(argc,argv, "-input", StringParm, infile, 1, use)
    CommandLineParm(argc,argv, "-format", StringParm, format, 1, use)
    CommandLineParm(argc,argv, "-output", StringParm, outfile, 1,use)
    CommandLineParm(argc,argv, "-xsize", IntegerParm, &xsize, 1, use)
    CommandLineParm(argc,argv, "-ysize", IntegerParm, &ysize, 1, use)
    CommandLineParm(argc, argv, "-iter", IntegerParm, &num_iter, 
                                    0, use)
    if (num_iter < 0) num_iter = 1
    CommandLineParm(argc, argv, "-cycles", IntegerParm, &num_cycles,
                                    0, use)
    if (num_cycles < 0) num_cycles = 1

    if (Keyword(format, "complex8"))    in_format = 0
    else if (Keyword(format, "complex4"))    in_format = 1
    else if (Keyword(format, "byte"))    in_format = 2
    else if (Keyword(format, "float"))    in_format = 3
    else {
        fprintf(stderr, "Unrecognized format: %s", format)
        exit(BAD_PARAMETER)


    print("Input file =    %s", infile)
    print("Input file type = %s", format)
    print("Output file =    %s", outfile)
    print("File dimensions = %dx%d (cols x rows).", xsize, ysize)

    #    OPEN FILES, ALLOCATE MEMORY
    OpenFile(&ifp, infile, "r")
    OpenFile(&ofp, outfile, "w")
    AllocateFloat(&phase, xsize*ysize, "phase data")
    AllocateFloat(&soln, xsize*ysize, "solution array")

    #    READ AND PROCESS DATA
    print("Reading phase data...")
    GetPhase(in_format, ifp, infile, phase, xsize, ysize)

    #    UNWRAP
    print("Unwrapping...")
    AllocateFloat(&dx, xsize*ysize, "dx derivs")
    AllocateFloat(&dy, ysize*ysize, "dy derivs")
    DxPhaseGradient(phase, dx, xsize, ysize)
    DyPhaseGradient(phase, dy, xsize, ysize)
    for (k=0; k<xsize*ysize; k += 1)    soln[k] = 0.0
    UnweightedMultigridUnwrap(soln, dx, dy, xsize, ysize,
                                                        num_cycles, num_iter)
    print("\nFinished")
    PrintMinAndMax(xsize, ysize, soln, "solution")

    #    SAVE RESULT
    for (k=0; k<xsize*ysize; k += 1) 
        soln[k] *= TWOPI
    print("Saving unwrapped surface to file '%s'", outfile)
    WriteFloat(ofp, soln, xsize*ysize, outfile)
    delete(dx)
    delete(dy)
    delete(soln)

#
 * mainunwt.c -- unweighted least-squares phase unwrapping
 *                             by means of direct transforms
 *
 * Source code files required:
 *                    dct.c         extract.c                grad.c         laplace.c
 *         mainunwt.c         solncos.c                util.c

#include <stdio.h>
#include <math.h>
#include "solncos.h"
#include "laplace.h"
#include "file.h"
#include "pi.h"
#include "util.h"
#include "extract.h"

main (int argc, char *argv[])
{
    int                        i, j, k, n
    FILE                     *ifp, *ofp
    float                    *phase;         # array
    float                    *soln;            # array
    double                 *xcos, *ycos
    char                     infile[200], outfile[200]
    char                     format[200]
    int                        in_format
    int                        xsize, ysize;     # dimensions of arrays
    int                        xsize_dct, ysize_dct
    double                 rmin, rmax, rscale
    double                 one_over_twopi = 1.0/TWOPI
    char                     use[] =             # define usage statement
     "Usage: prog-name -input file -format fkey -output file"
     "    -xsize x -ysize y"
     "where 'fkey' is a keyword designating the input file type"
     "(key = complex8, complex4, float or byte).    The dimensions"
     "x and y must be powers of two plus 1 (e.g., 257, 513, etc.)"

    print("Unweighted Phase Unwrapping by DCT/FFT")

    # GET COMMAND LINE PARAMETERS AND CHECK
    CommandLineParm(argc,argv, "-input", StringParm, infile, 1, use)
    CommandLineParm(argc,argv, "-format", StringParm, format, 1, use)
    CommandLineParm(argc,argv, "-output", StringParm, outfile, 1,use)
    CommandLineParm(argc,argv, "-xsize", IntegerParm, &xsize, 1, use)
    CommandLineParm(argc,argv, "-ysize", IntegerParm, &ysize, 1, use)

    if (Keyword(format, "complex8"))    in_format = 0
    else if (Keyword(format, "complex4"))    in_format = 1
    else if (Keyword(format, "byte"))    in_format = 2
    else if (Keyword(format, "float"))    in_format = 3
    else {
        fprintf(stderr, "Unrecognized format: %s", format)
        exit(BAD_PARAMETER)


    print("Input file =    %s", infile)
    print("Input file type = %s", format)
    print("Output file =    %s", outfile)
    print("File dimensions = %dx%d (cols x rows).", xsize, ysize)

    # Verify dimensions are of form 2**n + 1
    for (xsize_dct = 1; xsize_dct + 1 < xsize; xsize_dct*=2)
        
    xsize_dct = xsize_dct + 1
    for (ysize_dct = 1; ysize_dct + 1 < ysize; ysize_dct*=2)
        
    ysize_dct = ysize_dct + 1
    if (xsize_dct != xsize or ysize_dct != ysize) {
        fprintf(stderr, "Error: dims %dx%d must be 2**n + 1",
                        xsize, ysize)
        exit(BAD_PARAMETER)


    #    OPEN FILES, ALLOCATE MEMORY
    OpenFile(&ifp, infile, "r")
    OpenFile(&ofp, outfile, "w")
    AllocateFloat(&phase, xsize*ysize, "phase data")
    AllocateFloat(&soln, xsize*ysize, "scratch data")

    #    READ AND PROCESS DATA
    print("Reading phase data...")
    GetPhase(in_format, ifp, infile, phase, xsize, ysize)

    #    UNWRAP
    print("Computing Laplacian...")
    ComputeLaplacian(phase, soln, NULL, NULL, xsize, ysize, 1)
    delete(phase);    
    AllocateDouble(&xcos, xsize, "cosine terms")
    AllocateDouble(&ycos, ysize, "cosine terms")
    print("Performing direct transform...")
    DirectSolnByCosineTransform(soln, xsize, ysize, xcos, ycos)
    delete(xcos)
    delete(ycos)
    print("\nFinished")
    PrintMinAndMax(xsize, ysize, soln, "solution")

    #    SAVE RESULT
    for (k=0; k<xsize*ysize; k += 1)
        soln[k] *= TWOPI
    print("Saving unwrapped surface to file '%s'", outfile)
    WriteFloat(ofp, soln, xsize*ysize, outfile)
    delete(soln)

#
 *     maskcut.c -- function for quality map guided mask cut generation

#include <stdio.h>
#include <math.h>
#include "list.h"
#include "util.h"
#include "maskcut.h"

# Main function for Mask Cut Algorithm for phase unwrapping
int QualityGuidedMask(float *phase, float *qual_map, 
                                            unsigned char *bitflags, int xsize, int ysize, 
                                            int mask_code, UnwrapMode unwrap_mode,
                                            int debug_flag, char *infile)
{
    int    i, j, k, a, b, c, n, num_pieces=0
    float    value
    float    min_qual, small_val = -1.0E+10
    int        num_index, max_list_size, bench, benchout
    int        *index_list
    int        postponed_code=POSTPONED
    int        avoid_code
    int        charge, num_pixels

    bench = xsize*ysize/100
    if (bench < 1) bench = 1
    benchout = 10*bench
    max_list_size = xsize*ysize/2; # this size can be reduced
    AllocateInt(&index_list, max_list_size + 1,
                            "bookkeeping list (index)")

    # find starting point
    n = 0
    num_index = 0
    avoid_code = BORDER | mask_code
    for (j=0; j<ysize; j += 1) {
        for (i=0; i<xsize; i += 1) {
            min_qual = small_val
            num_index = 0
            k = j*xsize + i
            if ((bitflags[k] & RESIDUE) and !(bitflags[k] & avoid_code)) {
                charge = (bitflags[k] & POS_RES) ? 1 : -1
                 += 1num_pieces
                num_pixels = 1
                bitflags[k] |= mask_code
                if ((bitflags[k] & BORDER) 
                        or i==0 or i==xsize-1 or j==0 or j==ysize-1) {
                    continue

                UpdateList(qual_map, i, j, value, phase, NULL, bitflags,
                                xsize, ysize, index_list, &num_index,
                                postponed_code, VISITED, postponed_code,
                                max_list_size, unwrap_mode==dxdygrad, &min_qual)
                while (charge != 0 and num_index > 0) {
                    if (n%bench==0) {
                        print("%d ", n/bench)
                        fflush(stdout)

                     += 1n
                    if (!GetNextOneToUnwrap(&a, &b, index_list, &num_index,
                                                                    xsize, ysize))
                        break;     # no more to unwrap
                    c = b*xsize + a
                     += 1num_pixels; 
                    if ((bitflags[c] & RESIDUE) and !(bitflags[c] & mask_code)) 
                        charge += (bitflags[c] & POS_RES) ? 1 : -1
                    if ((bitflags[c] & BORDER) 
                            or a==0 or a==xsize-1 or b==0 or b==ysize-1) {
                        charge = 0

                    bitflags[c] |= mask_code
                    UpdateList(qual_map, a, b, value, phase, NULL, bitflags,
                                xsize, ysize, index_list, &num_index,
                                postponed_code, VISITED, postponed_code,
                                max_list_size, unwrap_mode==dxdygrad, &min_qual)

                # reset pixels on list
                for (k=0; k<num_index; k += 1) 
                    bitflags[index_list[k]] &= ~(VISITED)
                for (k=0; k<xsize*ysize; k += 1) 
                    bitflags[k] &= ~(VISITED)



    print("")
    delete(index_list)
    return num_pieces

#
 * maskfat.c -- fatten mask

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "util.h"
"""

def fattenmask(array, mask_code, thickness, xsize, ysize):
    """
    Fattens mask and border by "thickness" pixels
    """
    t2 = thickness
    if t2 < 1:
        return
    temp = xsize * ysize
    for j in xrange(ysize):
        for i in (xsize):
            k = j * xsize + i
            temp[k] = 0
            jj = j - t2
            # FIXME
            while True:
                if not(jj <= j + t2 and not temp[k]):
                    break
                jj += 1
                ii = i - t2
                while True:
                    if not(ii <= i + t2 and not(temp[k])):
                        break
                    ii += 1
                    kk = jj * xsize + ii
                    if ii >= 0 and ii < xsize and jj >= 0 and jj < ysize and (array[kk] & mask_code):
                        temp[k] = 1
    for k in xrange(xsize * ysize):
        if temp[k]:
            array[k] |= mask_code
    delete(temp)

"""

# Fattens zero values in quality map by "thickness" pixels
void FattenQual(float *qual_map, int thickness,
                                int xsize, int ysize)
{
    int i, j, k, ii, jj, kk, t2=thickness
    unsigned char *temp
    if (t2 < 1) return
    AllocateByte(&temp, xsize*ysize, "temp array")
    for (j=0; j<ysize; j += 1) {
        for (i=0; i<xsize; i += 1) {
            k = j*xsize + i
            temp[k] = 0
            for (jj=j-t2; jj<=j+t2 and !temp[k]; jj += 1) {
                for (ii=i-t2; ii<=i+t2 and !temp[k]; ii += 1) {
                    kk = jj*xsize + ii
                    if (ii>=0 and ii<xsize and jj>=0 and jj<ysize 
                                            and qual_map[kk]==0.0)
                        temp[k] = 1




    for (k=0; k<xsize*ysize; k += 1) {
        if (temp[k]) qual_map[k] = 0.0

    delete(temp)

#
 *     maskthin.c -- function for thinning a mask

#include <stdio.h>
#include <math.h>
#include "list.h"
#include "maskthin.h"

# Thins mask without disconnecting it.
void ThinMask(unsigned char *bitflags, int xsize, int ysize,
                            int cut_code)
{
    int     i, j, k, n, ii, jj, kk, k1, k2, k3
    int     numflips, removed_some=1, avoid_code
    avoid_code = RESIDUE | BORDER
    while (removed_some) { 
        removed_some = 0
        for (j=1; j<ysize - 1; j += 1) {
            for (i=1; i<xsize - 1; i += 1) {
                k = j*xsize + i
                # only consider branch cut pixels for removal
                if (!(bitflags[k] & cut_code)) continue
                # don't remove pixel k if it borders a residue or mask
                if (bitflags[k-xsize-1] & avoid_code) continue
                if (bitflags[k-xsize] & avoid_code) continue
                if (bitflags[k-xsize+1] & avoid_code) continue
                if (bitflags[k-1] & avoid_code) continue
                if (bitflags[k] & avoid_code) continue
                if (bitflags[k+1] & avoid_code) continue
                if (bitflags[k+xsize-1] & avoid_code) continue
                if (bitflags[k+xsize] & avoid_code) continue
                if (bitflags[k+xsize+1] & avoid_code) continue
                # walk around border of 3x3 nbhd; make sure removal
                # does not leave branch cut pixels disconnected in
                # 3x3 neighborhood
                kk = k - xsize;    # start at upper middle pixel
                for (n=0, numflips=0; n<4; n += 1) {
                    if (n==0) { k1 = k-xsize; k2 = k1+1; k3 = k2+xsize; }    
                    else if (n==1) { k1 = k+1; k2 = k1+xsize; k3 = k2-1; }    
                    else if (n==2) { k1 = k+xsize; k2 = k1-1; k3 = k2-xsize; }    
                    else if (n==3) { k1 = k-1; k2 = k1-xsize; k3 = k2+1; }    
                    if (bitflags[k1] & cut_code) { 
                        if (!(bitflags[k3] & cut_code))    
                             += 1numflips

                    else {
                        if (bitflags[k2] & cut_code) {
                             += 1numflips
                            if (!(bitflags[k3] & cut_code))  += 1numflips

                        else if (bitflags[k3] & cut_code)  += 1numflips


                if (numflips < 3) { # remove pixel from branch cut
                    if (!(bitflags[k-xsize]&cut_code) 
                                 or !(bitflags[k+1]&cut_code) 
                                            or !(bitflags[k+xsize]&cut_code) 
                                                     or !(bitflags[k-1]&cut_code)) {
                        bitflags[k] &= (~cut_code)
                         += 1removed_some




        print("Trimmed %d mask pixels", removed_some)

"""

"""
                      FUNCTION FOR UNWRAPPING AROUND CUTS
                      ===================================
"""

def unwraparoundcuts(phase, bitflags, soln, xsize, ysize, cut_code, debug_flag,
    infile):
    """
    Unwrap the phase data (by Itoh's method) without crossing any branch cuts.
    Return number of disconnected pieces.
    """
    num_pieces = 0
    small_val = -1.0E+10
    unwrapped_code = "UNWRAPPED"
    postponed_code = "POSTPONED"
    bench = xsize * ysize / 100
    if bench < 1:
        bench = 1
    benchout = 10 * bench
    min_qual = small_val
    max_list_size = xsize * ysize; # this size may be reduced
    index_list = max_list_size + 1 # bookkeeping list (index)
    avoid_code = cut_code | unwrapped_code | BORDER

    # find starting point
    n = 0
    num_index = 0
    for j in xrange(ysize):
        for i in xrange(xsize):
            k = j * xsize + i
            if not(bitflags[k] & avoid_code):
                bitflags[k] |= unwrapped_code
                if bitflags[k] & postponed_code:
                    # soln[k] already stores the unwrapped value
                    value = soln[k]
                else:
                    num_pieces += 1
                    value = soln[k] = phase[k]
                updatelist(NULL, i, j, value, phase, soln, bitflags, xsize,
                    ysize, index_list, num_index, avoid_code, unwrapped_code,
                    postponed_code, max_list_size, 0, min_qual)
                while num_index > 0:
                    if n % bench == 0:
                        print("%d ", n/bench)
                        if 0 and debug_flag and n % benchout == 0 and n > 0:
                            file_num = 0
                            filename = "%s.fill-%d" % (infile, file_num)
                            savebytetoimage(bitflags, "fill path", filename,
                                xsize, ysize, 1, 1, avoid_code)
                    n += 1
                    if not(getnextonetounwrap(a, b, index_list, num_index, xsize,
                        ysize)):
                        break # no more to unwrap
                    c = b * xsize + a
                    bitflags[c] |= unwrapped_code
                    value = soln[c]                
                    updatelist(a, b, value, phase, soln, bitflags, xsize, ysize,
                        index_list, num_index, avoid_code, unwrapped_code,
                        postponed_code, max_list_size, 0, min_qual)
    print("")
    delete(index_list)
    # unwrap branch cut pixels
    for j in xrange(1, ysize):
        for i in xrange(1, xsize):
            k = j * xsize + i
            if bitflags[k] & cut_code:
                if not(bitflags[k-1] & cut_code):
                    soln[k] = soln[k - 1] + gradient(phase[k], phase[k-1])
                elif not(bitflags[k - xsize] & cut_code):
                    soln[k] = soln[k-xsize] + gradient(phase[k], phase[k-xsize])
    return num_pieces

"""
#
 *    pcg.c -- function for Preconditioned Conjugate Gradient solution
 *                     of weighted least-squares phase-unwrapping problem

#include <stdio.h>
#include <math.h>
#include "solncos.h"
#include "pcg.h"
#include "util.h"
#define SIMIN(x,y) (((x)*(x) < (y)*(y)) ? (x)*(x) : (y)*(y))

# Main function of PCG algorithm for phase unwrapping.    If dywts
# is null, then the dx and dy weights are calculated from dxwts
void PCGUnwrap(float *rarray, float *zarray, float *parray,
                             float *soln, float *dxwts, float *dywts, int xsize,
                             int ysize, int max_iter, double epsi_con)
{
    int            i, j, k, iloop
    double     sum, alpha, beta, beta_prev, epsi
    double     *xcos, *ycos
    for (k=0, sum=0.0; k<xsize*ysize; k += 1)
        sum += rarray[k]*rarray[k]
    sum = sqrt(sum/(xsize*ysize))
    AllocateDouble(&xcos, xsize, "cosine terms")
    AllocateDouble(&ycos, ysize, "cosine terms")
    for (iloop=0; iloop < max_iter; iloop += 1) {
        PCGIterate(rarray, zarray, parray, soln, dxwts, dywts, 
                             xsize, ysize, xcos, ycos, iloop, sum, &alpha,
                             &beta, &beta_prev, &epsi)
        if (epsi < epsi_con) {
            print("Breaking out of main loop (due to convergence)")
            break


    delete(xcos)
    delete(ycos)
 

# Main function within iterative loop of PCG algorithm.
void PCGIterate(float *rarray, float *zarray, float *parray,
                                float *soln, float *dxwts, float *dywts, int xsize,
                                int ysize, double *xcos, double *ycos, int iloop,
                                double sum0, double *alpha, double *beta,
                                double *beta_prev, double *epsi)
{
    int         i, j, k
    int         k1, k2, k3, k4
    double    sum, w1, w2, w3, w4, btemp, delta, avg, scale
    float     *wts

    scale = 1.0/(xsize*ysize)
    # remove constant bias from rarray
    for (k=0, avg=0.0; k<xsize*ysize; k += 1)    avg += rarray[k]
    avg *= scale
    for (k=0; k<xsize*ysize; k += 1)    rarray[k] -= avg
    # compute cosine transform solution of Laplacian in rarray
    for (k=0; k<xsize*ysize; k += 1)    {
        zarray[k] = rarray[k]

    DirectSolnByCosineTransform(zarray, xsize, ysize, xcos, ycos)
    # calculate beta and parray
    for (k=0, *beta=0.0; k<xsize*ysize; k += 1) {
        *beta += rarray[k]*zarray[k]

    print("beta = %lf", *beta)
    if (iloop == 0) {
        for (k=0; k<xsize*ysize; k += 1) {
            parray[k] = zarray[k]


    else {
        btemp = (*beta)/(*beta_prev)
        for (k=0; k<xsize*ysize; k += 1) {
            parray[k] = zarray[k] + btemp*parray[k]


    # remove constant bias from parray
    for (k=0, avg=0.0; k<xsize*ysize; k += 1)    avg += parray[k]
    avg *= scale
    for (k=0; k<xsize*ysize; k += 1)    parray[k] -= avg
    *beta_prev = *beta
    # calculate Qp
    for (j=0; j<ysize; j += 1) {
        for (i=0; i<xsize; i += 1) {
            k = j*xsize + i
            k1 = (i<xsize-1) ? k + 1 : k - 1
            k2 = (i>0) ? k - 1 : k + 1
            k3 = (j<ysize-1) ? k + xsize : k - xsize
            k4 = (j>0) ? k - xsize : k + xsize
            if (dxwts==NULL and dywts==NULL) {    # unweighted
                w1 = w2 = w3 = w4 = 1.0

            else if (dxwts==NULL or dywts==NULL) {    # one set of wts
                wts = (dxwts) ? dxwts : dywts
                w1 = SIMIN(wts[k], wts[k1])
                w2 = SIMIN(wts[k], wts[k2])
                w3 = SIMIN(wts[k], wts[k3])
                w4 = SIMIN(wts[k], wts[k4])

            else {        # dxwts and dywts are both supplied
                w1 = dxwts[k]
                w2 = (i>0) ? dxwts[k-1] : dxwts[k]
                w3 = dywts[k]
                w4 = (j>0) ? dywts[k-xsize] : dywts[k]


            zarray[k] = (w1 + w2 + w3 + w4)*parray[k]
                                                - (w1*parray[k1] + w2*parray[k2] 
                                                                + w3*parray[k3] + w4*parray[k4])


    # calculate alpha
    for (k=0, *alpha=0.0; k<xsize*ysize; k += 1) {
        *alpha += zarray[k]*parray[k]

    *alpha = *beta/(*alpha)
    print("alpha = %lf", *alpha)
    # update rarray
    for (k=0; k<xsize*ysize; k += 1) {
        rarray[k] -= (*alpha)*zarray[k]

    # update parray
    for (k=0; k<xsize*ysize; k += 1) {
        soln[k] += (*alpha)*parray[k]

    # remove constant bias from soln
    for (k=0, avg=0.0; k<xsize*ysize; k += 1)    avg += soln[k]
    avg *= scale
    for (k=0; k<xsize*ysize; k += 1)    soln[k] -= avg
    # compute epsi and delta
    for (k=0, sum=0.0; k<xsize*ysize; k += 1) {
        sum += rarray[k]*rarray[k]

    *epsi = sqrt(sum/(xsize*ysize))/sum0
    for (k=0, sum=0.0; k<xsize*ysize; k += 1) {
        sum += (*alpha)*(*alpha)*parray[k]*parray[k]

    delta = sqrt(sum/(xsize*ysize))
    print("ITER %d: EPSI %lf DELTA %lf", iloop, *epsi, delta)

#
 *    qualgrad.c -- functions for computing max gradient quality map

#include <stdio.h>
#include <math.h>
#include "qualgrad.h"
#include "dxdygrad.h"

# Compute max gradient quality map in array "result".    Size
# of averaging template is tsize x tsize pizels.    The array
# temp1 is required for temporary data.
void MaxPhaseGradients(float *phase, float *result,
                                        unsigned char *bitflags, int ignore_code,
                                        float *temp1, int tsize, int xsize, int ysize)
{
    int    add_flag
    print("Extracting dx gradients")
    DxPhaseGradient(phase, temp1, xsize, ysize)
    add_flag = 0
    print("Extracting dx max's")
    DxGradMax(temp1, result, xsize, ysize, tsize, bitflags,
                        ignore_code, add_flag)
    print("Extracting dy gradients")
    DyPhaseGradient(phase, temp1, xsize, ysize)
    add_flag = 1
    print("Extracting dy max's")
    # add to dx max's
    DxGradMax(temp1, result, xsize, ysize, tsize, bitflags,
                        ignore_code, add_flag)


# Compute max dx gradients in array "dxmax". If add_code
# is 1, the results are added to the dxmax array,
# otherwise they overwrite.
void DxGradMax(float *dx, float *dxmax, int xsize, int ysize,
     int tsize, unsigned char *bitflags, int avoid_code, int add_code)
{
    int        i, j, a, b, c, hs
    float    r, dmax
    if (tsize < 3) {
        for (i=0; i<xsize*ysize; i += 1)
            dxmax[i] = (add_code) ? dxmax[i] + dx[i] : dx[i]

    else {
        hs = tsize/2
        for (j=0; j<ysize; j += 1) {
            for (i=0; i<xsize; i += 1) {
                for (dmax = 0.0, b=j-hs; b<=j+hs; b += 1) {
                    if (b < 0 or b >= ysize) continue
                    for (a=i-hs; a < i+hs; a += 1) {    # a < to match dx
                        if (a < 0 or a >= xsize) continue
                        c = b*xsize + a
                        if (bitflags and (bitflags[c]&avoid_code)) continue
                        r = dx[c]
                        if (r < 0) r = -r
                        if (dmax < r) dmax = r; 


                if (add_code)
                    dxmax[j*xsize + i] += dmax
                else
                    dxmax[j*xsize + i] = dmax

     

 

# Compute max dy gradients in array "dymax". If add_code
# is 1, the results are added to the dymax array,
# otherwise they overwrite.
void DyGradMax(float *dy, float *dymax, int xsize, int ysize,
     int tsize, unsigned char *bitflags, int avoid_code, int add_code)
{
    int        i, j, a, b, c, hs
    float    r, dmax
    if (tsize < 3) {
        for (i=0; i<xsize*ysize; i += 1)
            dymax[i] = (add_code) ? dymax[i] + dy[i] : dy[i]

    else {
        hs = tsize/2
        for (j=0; j<ysize; j += 1) {
            for (i=0; i<xsize; i += 1) {
                for (dmax = 0.0, b=j-hs; b < j+hs; b += 1) {    
                    # b < j+hs to match dy
                    if (b < 0 or b >= ysize) continue
                    for (a=i-hs; a<=i+hs; a += 1) {
                        if (a < 0 or a >= xsize) continue
                        c = b*xsize + a
                        if (bitflags and (bitflags[c] & avoid_code)) continue
                        r = dy[c]
                        if (r < 0) r = -r
                        if (dmax < r) dmax = r; 


                if (add_code)
                    dymax[j*xsize + i] += dmax
                else
                    dymax[j*xsize + i] = dmax

     

 
#
 *    quality.c -- function for quality-guided path-following

#include <stdio.h>
#include <math.h>
#include "quality.h"
#include "list.h"
#include "util.h"

# Main function for quality-guided algorithm for phase unwrapping
int QualityGuidedPathFollower(float *phase, float *qual_map,
                                                        unsigned char *bitflags, float *soln, 
                                                        int xsize, int ysize, int avoid_code, 
                                                        UnwrapMode unwrap_mode, int debug_flag,
                                                        char *infile)
{
    int    i, j, k, a, b, c, n, num_pieces=0
    float    value
    float    min_qual, small_val = -1.0E+10
    int        num_index, max_list_size, bench, benchout
    int        *index_list
    int        postponed_code=POSTPONED, unwrapped_code=UNWRAPPED

    bench = xsize*ysize/100
    if (bench < 1) bench = 1
    benchout = 10*bench
    min_qual = small_val
    max_list_size = (xsize + ysize)
    AllocateInt(&index_list, max_list_size + 1, 
                            "bookkeeping list (index)")
    avoid_code |= unwrapped_code
    # find starting point
    n = 0
    num_index = 0
    for (j=0; j<ysize; j += 1) {
        for (i=0; i<xsize; i += 1) {
            k = j*xsize + i
            if (!(bitflags[k] & avoid_code)) {
                if (bitflags[k] & postponed_code) 
                    # soln[k] already stores the unwrapped value
                    value = soln[k]
                else {
                     += 1num_pieces
                    value = soln[k] = phase[k]
 
                bitflags[k] |= unwrapped_code
                bitflags[k] &= ~postponed_code
                UpdateList(qual_map, i, j, value, phase, soln, bitflags,
                                     xsize, ysize, index_list, &num_index, avoid_code, 
                                     unwrapped_code, postponed_code, max_list_size,
                                     unwrap_mode==dxdygrad, &min_qual)
                while (num_index > 0) {
                    if (n%bench==0) {
                        print("%d ", n/bench)
                        fflush(stdout)
                        if (0 and debug_flag and n%benchout==0 and n>0) {
                            char filename[300]
                            static char file_num = 0
                            sprintf(filename, "%s.path-%d", infile,  += 1file_num)
                            SaveByteToImage(bitflags, "unwrapping path", filename,
                                                            xsize, ysize, 1, 1, avoid_code)


                     += 1n
                    if (!GetNextOneToUnwrap(&a, &b, index_list, &num_index,
                                xsize, ysize)) break;     # no more to unwrap
                    c = b*xsize + a
                    bitflags[c] |= unwrapped_code
                    bitflags[c] &= ~postponed_code
                    value = soln[c]
                    UpdateList(qual_map, a, b, value, phase, soln, bitflags,
                                     xsize, ysize, index_list, &num_index, avoid_code, 
                                     unwrapped_code, postponed_code, max_list_size,
                                     unwrap_mode==dxdygrad, &min_qual)
                    if (num_index <= 0) {
                        min_qual = small_val
                        for (c=0; c<xsize*ysize; c += 1) {
                            if ((bitflags[c] & postponed_code) 
                                             and !(bitflags[c] & unwrapped_code)) {
                                a = c%xsize;        
                                b = c/xsize
                                InsertList(soln, soln[c], qual_map, bitflags, a, b,
                                                     index_list, &num_index, xsize, ysize,
                                                     unwrapped_code, postponed_code,
                                                     &min_qual, max_list_size)


                        if (num_index > 0) min_qual = qual_map[index_list[0]]

    # while ...

     # for (b ...
     # for (a ...
    print("")
    delete(index_list)
    return num_pieces

#
 *     qualpseu.c -- functions for computing pseudocorrelation 
 *                                 quality map

#include <stdio.h>
#include <math.h>
#include "pi.h"
#include "qualpseu.h"
#include "dxdygrad.h"

# Compute pseudocorrelation quality map in array "result".
# The size of the averaging template is tsize x tsize pizels.
# The array temp1 is required for temporary data.
void PseudoCorrelation(float *phase, float *result,
                                         unsigned char *bitflags, int ignore_code,
                                         float *temp1, int tsize, int xsize, int ysize)
{
    int    i, j, k, add_flag
    print("Extracting cos's")
    for (k=0; k<xsize*ysize; k += 1) 
        temp1[k] = cos(TWOPI*phase[k])
    add_flag = 0
    print("Filtering cos's")
    SqrAvgFilter(temp1, result, xsize, ysize, tsize,
                             bitflags, ignore_code, add_flag)
    print("Extracting sin's")
    for (k=0; k<xsize*ysize; k += 1) 
        temp1[k] = sin(TWOPI*phase[k])
    add_flag = 1
    print("Filtering sin's")
    SqrAvgFilter(temp1, result, xsize, ysize, tsize,
                             bitflags, ignore_code, add_flag)
    print("Square root")
    for (k=0; k<xsize*ysize; k += 1) 
        result[k] = 1.0 - sqrt(result[k])


# Computes average of array 'in' in tsize x tsize windows,
# then squares the values.    If add_code is 1, the results
# are added to the values in out, otherwise they overwrite.
void SqrAvgFilter(float *in, float *out, int xsize, int ysize,
     int tsize, unsigned char *bitflags, int avoid_code, int add_code)
{
    int        i, j, k, a, b, aa, bb, cc, n, hs
    float    r, avg
    if (tsize < 3 and !add_code) {
        for (i=0; i<xsize*ysize; i += 1)
            out[i] = 0.0

    else {
        hs = tsize/2
        for (j=0; j<ysize; j += 1) {
            for (i=0; i<xsize; i += 1) {
                avg = 0.0
                for (n=0, b=j-hs; b<=j+hs; b += 1) {
                    if ((bb = b) < 0) bb = -bb
                    else if (bb >= ysize) bb = 2*ysize - 2 - bb; 
                    for (a=i-hs; a <= i+hs; a += 1) {    
                        if ((aa = a) < 0) aa = -aa
                        else if (aa >= xsize) aa = 2*xsize - 2 - aa; 
                        cc = bb*xsize + aa
                        if (aa>=0 and aa<xsize-1 and bb>=0 and bb<ysize-1) {
                            r = in[cc]
                            avg += r
                             += 1n;    



                r = (n>0) ? 1.0/n : 0.0
                avg *= r
                if (add_code)
                    out[j*xsize + i] += avg*avg
                else
                    out[j*xsize + i] = avg*avg; 

     

 
#
 *    qualvar.c -- functions for computing phase derivative variance
 *                             quality map

#include <stdio.h>
#include <math.h>
#include "qualvar.h"
#include "dxdygrad.h"

# Compute variance of phase derivatives in array "result".
# The size of the averaging template is tsize x tsize pizels.
# The array temp1 is required for temporary data.
void PhaseVariance(float *phase, float *result, 
                                     unsigned char *bitflags, int ignore_code,
                                     float *temp1, int tsize, int xsize, int ysize)
{
    int    add_flag
    print("Extracting dx gradients")
    DxPhaseGradient(phase, temp1, xsize, ysize)
    add_flag = 0
    print("Extracting dx variances")
    DxGradVar(temp1, result, xsize, ysize, tsize, bitflags,
                        ignore_code, add_flag)
    print("Extracting dy gradients")
    DyPhaseGradient(phase, temp1, xsize, ysize)
    add_flag = 1
    print("Extracting dy variances")
    # add to dx variances
    DyGradVar(temp1, result, xsize, ysize, tsize, bitflags,
                        ignore_code, add_flag)


# Compute variance of dx phase derivs in array "dxvar" in
# tsize x tsize windows.    If add_code is 1, the results
# are added to the dxvar array, otherwise they overwrite.
void DxGradVar(float *dx, float *dxvar, int xsize, int ysize,
    int tsize, unsigned char *bitflags, int avoid_code, int add_code)
{
    int        i, j, k, a, b, aa, bb, cc, n, hs
    float    r, avg, avgsqr
    if (tsize < 3 and !add_code) {
        for (i=0; i<xsize*ysize; i += 1)
            dxvar[i] = 0.0

    else {
        hs = tsize/2
        for (j=0; j<ysize; j += 1) {
            for (i=0; i<xsize; i += 1) {
                avg = avgsqr = 0.0
                for (n=0, b=j-hs; b<=j+hs; b += 1) {
                    if ((bb = b) < 0) bb = -bb
                    else if (bb >= ysize) bb = 2*ysize - 2 - bb; 
                    for (a=i-hs; a < i+hs; a += 1) {     # a < to match dx
                        if ((aa = a) < 0) aa = -aa
                        else if (aa >= xsize) aa = 2*xsize - 2 - aa; 
                        cc = bb*xsize + aa
                        if (bitflags and (bitflags[cc]&avoid_code)) continue
                        r = dx[cc]
                        avg += r
                        avgsqr += r*r
                         += 1n;    


                r = (n>0) ? 1.0/n : 0.0
                avg *= r
                avgsqr *= r
                if (add_code)
                    dxvar[j*xsize + i] += avgsqr - avg*avg;     # variance
                else
                    dxvar[j*xsize + i] = avgsqr - avg*avg;     # variance

     

 

# Compute variance of dy phase derivs in array "dyvar" in
# tsize x tsize windows.    If add_code is 1, the results
# are added to the dyvar array, otherwise they overwrite.
void DyGradVar(float *dy, float *dyvar, int xsize, int ysize,
    int tsize, unsigned char *bitflags, int avoid_code, int add_code)
{
    int        i, j, k, a, b, aa, bb, cc, n, hs
    float    r, avg, avgsqr
    if (tsize < 3 and !add_code) {
        for (i=0; i<xsize*ysize; i += 1)
            dyvar[i] = 0.0

    else {
        hs = tsize/2
        for (j=0; j<ysize; j += 1) {
            for (i=0; i<xsize; i += 1) {
                avg = avgsqr = 0.0
                for (n=0, b=j-hs; b < j+hs; b += 1) { # b < to match dy
                    if ((bb = b) < 0) bb = -bb
                    else if (bb >= ysize) bb = 2*ysize - 2 - bb; 
                    for (a=i-hs; a<=i+hs; a += 1) {
                        if ((aa = a) < 0) aa = -aa
                        else if (aa >= xsize) aa = 2*xsize - 2 - aa; 
                        cc = bb*xsize + aa
                        if (bitflags and (bitflags[cc]&avoid_code)) continue
                        r = dy[cc]
                        avg += r
                        avgsqr += r*r
                         += 1n


                r = (n>0) ? 1.0/n : 0.0
                avg *= r
                avgsqr *= r
                if (add_code)
                    dyvar[j*xsize + i] += avgsqr - avg*avg;     # variance
                else
                    dyvar[j*xsize + i] = avgsqr - avg*avg;     # variance

     

 
#
 *    raster.c -- simple raster function for unwrapping residue-free
 *                            phase data by means of Itoh's method.

#include <stdio.h>
#include <math.h>
#include "raster.h"
#include "grad.h"
# Unwrap phase using simple raster technique based on Itoh's
# method.    This method works inly if there are no residues.
void RasterUnwrap(float *phase, float *soln, int xsize, int ysize)
{
    int        i, j, k
    soln[0] = phase[0]
    for (j=1; j<ysize; j += 1) {    # unwrap first column
        k = j*xsize + 0
        soln[k] = soln[k-xsize] + Gradient(phase[k], phase[k-xsize])

    for (j=0; j<ysize; j += 1) {     # unwrap all the rows
        for (i=1; i<xsize; i += 1) {
            k = j*xsize + i
            soln[k] = soln[k-1] + Gradient(phase[k], phase[k-1])



# 
 * relax.c - function for solving weighted or unweighted
 *                     version of Poisson equation (i.e., weighted
 *                     or unweighted least-squares phase unwrapping)
 *                     by means of red-black Gauss-Seidel relaxation

#include <stdio.h>
#include <math.h>
#include "relax.h"

#
 * Black-red Gauss-Seidel relaxation: assumes rho(i,j) 
 *                = 0.25*(dx(i,j) - dx(i+1,j) + dy(i,j) - dy(i,j+1)
 * Thus, f(i,j) 
 *        = 0.25*(f(i-1,j) + f(i+1,j) + f(i,j-1) + f(i,j+1)) + rho(i,j)
 * This version solves the weighted equation.    For the unweighted
 * equation, just pass null pointers for dxwts and dywts.

void Relax(float *soln, float *dx, float *dy, float *dxwts,
                     float *dywts, int w, int h, int numit)
{
    int i, j, k, n, ipass, istart
    float    x1, x2, x3, x4, y, z, r
    float    w1, w2, w3, w4, norm=0.0
    double    diff, maxdiff, sum

    if (w*h < 8*8 and numit < 2) numit = 2*(w + h)
    for (n=0; n<numit; n += 1) {                         # iteration
        for (ipass=0; ipass<2; ipass += 1) {     # red and black sweeps
            istart = ipass
            sum = maxdiff = 0.0
            for (j=0; j<h; j += 1) {                         # row
                for (i=istart; i<w; i+=2) {         # column
                    k = j*w + i
                    # weights
                    if (dxwts and dywts) {
                        w1 = (i < w - 1) ? dxwts[k+1] : dxwts[k-1]
                        if (w1 > dxwts[k]) w1 = dxwts[k]
                        w1 = w1*w1
                        w2 = (i > 0) ? dxwts[k-1] : dxwts[k+1];    
                        if (w2 > dxwts[k]) w2 = dxwts[k]
                        w2 = w2*w2
                        w3 = (j < h - 1) ? dywts[k+w] : dywts[k-w]
                        if (w3 > dywts[k]) w3 = dywts[k]
                        w3 = w3*w3
                        w4 = (j > 0) ? dywts[k-w] : dywts[k+w]
                        if (w4 > dywts[k]) w4 = dywts[k]
                        w4 = w4*w4
                        norm = w1 + w2 + w3 + w4
 
                    # surface points
                    x1 = (i < w - 1) ? soln[k+1] : soln[k-1]
                    x2 = (i > 0) ? soln[k-1] : soln[k+1]
                    x3 = (j < h - 1) ? soln[k+w] : soln[k-w]
                    x4 = (j > 0) ? soln[k-w] : soln[k+w]
                    # derivatives
                    y = (i > 0) ? dx[k-1] : -dx[k]
                    z = (j > 0) ? dy[k-w] : -dy[k]
                    if (norm > 1.0e-6) { 
                        r = (w1*(x1 - dx[k]) + w2*(x2 + y) 
                                                        + w3*(x3 - dy[k]) + w4*(x4 + z))/norm; 
 
                    else {
                        r = 0.25*(x1 - dx[k] + x2 + y + x3 - dy[k] + x4 + z)

                    diff = r - soln[k]
                    if (diff < 0) diff = -diff
                    maxdiff = (maxdiff < diff) ? diff : maxdiff
                    sum += diff*diff
                    soln[k] = r

                istart = (istart) ? 0 : 1



    for (n=1; n<w; n*=2) print("    ")
    print("Relax %d %d: rms = %lf", w, h, sqrt(sum/(w*h)))
"""

"""
                  FUNCTION FOR EXTRACTING RESIDUES FROM PHASE
                  ===========================================
"""

def residues(phase, bitflags, posres_code, negres_code, avoid_code, xsize, ysize):
    """
    Detect residues in phase data and mark them as positive or negative
    residues in the bitflags array. Ignore the pixels marked with avoid_code in
    the bitflags araay.
    """
    numres = 0
    for j in xrange(ysize - 1):
        for i in xrange(xsize - 1):
            k = j * xsize + i
            if (bitflags and (
                    (bitflags[k] & avoid_code) or (bitflags[k+1] & avoid_code)
                    or (bitflags[k + 1 + xsize] & avoid_code)
                    or (bitflags[k + xsize] & avoid_code))):
                continue; # masked region: don't unwrap
            r = (gradient(phase[k + 1], phase[k])
                + gradient(phase[k + 1 + xsize], phase[k + 1])
                + gradient(phase[k + xsize], phase[k + 1 + xsize])
                + gradient(phase[k], phase[k + xsize]))
            if bitflags:
                if r > 0.01:
                    bitflags[k] |= posres_code
                elif r < -0.01:
                    bitflags[k] |= negres_code
            if r ** 2 > 0.01:
                numres += 1
    return numres

"""
#
 *    solncos.c -- function for solving Poisson's equation by means of
 *                             fast cosine transform

#include <stdio.h>
#include <math.h>
#include "dct.h"
#include "pi.h"

# Main function for direct solution of Poisson's equation for
# unweighted least squares phase unwrapping by means of discrete
# cosine transform (which actually uses FFTs to compute the DCTs)
void DirectSolnByCosineTransform(float *array, int xsize, int ysize,
                                                                 double *xcos, double *ycos)
{
    int     i, j, m, n
    #    Transform
    print("Forward transform...")
    FastCosineTransform(array, xsize, ysize, 1)

    #    Divide
    print("Scaling...")
    for (i=0; i<xsize; i += 1) {
        xcos[i] = cos(i*PI/(xsize - 1.0))

    for (j=0; j<ysize; j += 1) {
        ycos[j] = cos(j*PI/(ysize - 1.0))

    array[0] = 0.0
    for (j=0; j<ysize; j += 1) {
        for (i=0; i<xsize; i += 1) {
            if (i==0 and j==0) {
                array[0] = 0.0

            else {
                array[j*xsize + i] /= (4.0 - 2.0*xcos[i] - 2.0*ycos[j])



    #    Transform
    print("Inverse transform...")
    FastCosineTransform(array, xsize, ysize, -1)

#
 *    trees.c -- functions for managing trees, loops, edges, etc.
 *                         in Flynn's min. discontinuity algorithm

#include <stdio.h>
#include <math.h>
#include "trees.h"

# Remove the loop and update all of the edges.
void RemoveLoop(int ibase, int jbase, int ilast, int jlast,
                                int *value, unsigned char *bitflags, short *vjump,
                                short *hjump, int xsize, int ysize)
{
    int     jtip,itip
    int     value_shift
    # Remove edge from last node to base of loop
    if (jbase>jlast)             --vjump[jbase*(xsize+1) + ibase]
    else if (jbase<jlast)     += 1vjump[jlast*(xsize+1) + ilast]
    else if (ibase>ilast)     += 1hjump[jbase*(xsize+1) + ibase]
    else if (ibase<ilast)    --hjump[jlast*(xsize+1) + ilast]
    # Back up over loop: remove edges and change subtree values
    jtip = jlast
    itip = ilast
    do {
        value_shift = -value[jtip*(xsize+1) + itip]
        ChangeOrphan(itip, jtip, value_shift, value, bitflags,
                                 xsize, ysize)
        if (jtip > 0 
                     and (bitflags[(jtip-1)*(xsize+1) + itip] & RIGHT)) {
            vjump[jtip*(xsize+1) + itip]--
            bitflags[(jtip-1)*(xsize+1) + itip] &= ~RIGHT
            jtip--

        else if ((jtip<ysize) &&
                (bitflags[(jtip+1)*(xsize+1) + itip] & LEFT)) {
            vjump[(jtip+1)*(xsize+1) + itip] += 1
            bitflags[(jtip+1)*(xsize+1) + itip] &= ~LEFT
            jtip += 1

        else if (itip > 0 
                                and (bitflags[jtip*(xsize+1) + itip-1] & DOWN)) {
            hjump[jtip*(xsize+1) + itip] += 1
            bitflags[jtip*(xsize+1) + itip-1] &= ~DOWN
            itip--

        else if ((itip<xsize) 
                                and (bitflags[jtip*(xsize+1) + itip+1] & UP)) {
            hjump[jtip*(xsize+1) + itip+1]--
            bitflags[jtip*(xsize+1) + itip+1] &= ~UP
            itip += 1

        else {
            print("Error: broken loop at col %d, row %d",itip,jtip)
            exit(-1)

 while ((jtip != jbase) or (itip != ibase))


#
 * Add 'value_change' (which is > 0) to the values of all nodes
 * of the tree rooted at the node (i,j).    Set a flag if the tree
 * includes (ilast,jlast). Make (i,j) child in the next iteration.

void ChangeExten(int i, int j, int ilast, int jlast,
                                 int *loop_found, int value_change, int *value,
                                 unsigned char *bitflags, int xsize, int ysize)
{
    unsigned char kids
    *loop_found |= ((j == jlast) and (i == ilast))
    kids = bitflags[j*(xsize+1) + i]
    if (kids & LEFT)    
        ChangeExten(i, j-1, ilast, jlast, loop_found, value_change,
                                value, bitflags, xsize, ysize)
    if (kids & RIGHT)    
        ChangeExten(i, j+1, ilast, jlast, loop_found, value_change,
                                value, bitflags, xsize, ysize)
    if (kids & UP)        
        ChangeExten(i-1, j, ilast, jlast, loop_found, value_change,
                                value, bitflags, xsize, ysize)
    if (kids & DOWN) 
        ChangeExten(i+1, j, ilast, jlast, loop_found, value_change,
                                value, bitflags, xsize, ysize)
    value[j*(xsize+1) + i] += value_change
    bitflags[j*(xsize+1) + i] |= (NEXT_TIME | THIS_TIME)


#
 * Add 'value_shift' (which is < 0) to the values of all nodes
 * of the tree rooted at the node (i,j).    If the change would make
 * 'value' < 0, make 'value' 0 instead    and split the node from
 * its parent.    Make (i,j) child in the next iteration.

void ChangeOrphan(int i, int j, int value_shift, int *value,
                                    unsigned char *bitflags, int xsize, int ysize)
{
    unsigned char kids
    if (value[j*(xsize+1) + i] + (value_shift) < 0) {
        value_shift = -value[j*(xsize+1) + i]
        if (j>0)         bitflags[(j-1)*(xsize+1) + i] &=    ~RIGHT
        if (j<ysize) bitflags[(j+1)*(xsize+1) + i] &=    ~LEFT
        if (i>0)         bitflags[j*(xsize+1) + i-1]     &=    ~DOWN
        if (i<xsize) bitflags[j*(xsize+1) + i+1]     &=    ~UP

    kids = bitflags[j*(xsize+1) + i]
    if (kids & LEFT) ChangeOrphan(i, j-1, value_shift, value,
                                                                bitflags, xsize, ysize)
    if (kids & RIGHT) ChangeOrphan(i, j+1, value_shift, value,
                                                                 bitflags, xsize, ysize)
    if (kids & UP) ChangeOrphan(i-1, j, value_shift, value,
                                                            bitflags, xsize, ysize)
    if (kids & DOWN) ChangeOrphan(i+1, j, value_shift, value,
                                                                bitflags, xsize, ysize)
    bitflags[j*(xsize+1) + i] |= (NEXT_TIME | THIS_TIME)
    value[j*(xsize+1) + i] += value_shift

#
 *    unfmg.c -- function for multigrid solution of UNWEIGHTED
 *                         least-squares phase-unwrapping problem

#include <stdio.h>
#include <math.h>
#include "unfmg.h"
#include "ungrid.h"
#include "dxdygrad.h"
# Unweighted multigrid function.
# Initialize soln array to zero before calling.
void UnweightedMultigridUnwrap(float *soln, float *dx, float *dy,
                                int xsize, int ysize, int num_cycles, int num_iter)
{
    int    n, coarsest_dim=3
    for (n=0; n<num_cycles; n += 1) {
        print("\nFMG CYCLE %d", n+1)
        Ufmg(soln, dx, dy, xsize, ysize, num_iter, coarsest_dim)


#
 * ungrid.c - functions for UNWEIGHTED multigrid phase
 *                        unwrapping 

#include <stdio.h>
#include <math.h>
#include "ungrid.h"
#include "gridops.h"
#include "gridmem.h"
#include "relax.h"

# Main function for unweighted multigrid phase unwrapping.
# (Called recursively.)
void Ufmg(float *soln, float *dx, float *dy,
                    int w, int h, int numit, int mindim)
{
    float    *dx2=NULL, *dy2=NULL, *soln2=NULL
    int        w2 = w/2, h2 = h/2
    if (!Coarsest(w, h, mindim)) {
        dx2 = Allocate(w2, h2, dx_type)
        dy2 = Allocate(w2, h2, dy_type)
        soln2 = Allocate(w2, h2, soln_type)
        Restrict(dx2, dy2, w2, h2, dx, dy, NULL, NULL, soln, w, h)
        Zero(soln2, w2, h2)
        Ufmg(soln2, dx2, dy2, w2, h2, numit, mindim)
        ProlongAndAccumulate(soln, w, h, soln2, w2, h2, NULL, NULL)

    # perform V-cycle multigrid on fine grid
    Umv(soln, dx, dy, w, h, numit, mindim)


# V-cycle function for unweighted multigrid phase unwrapping.
# (Called recursively.)
void Umv(float *soln, float *dx, float *dy,
                 int w, int h, int numit, int mindim)
{
    float *dx2=NULL, *dy2=NULL, *soln2=NULL
    int        w2 = w/2, h2 = h/2
    if (!Coarsest(w, h, mindim)) {
        Relax(soln, dx, dy, NULL, NULL, w, h, numit)
        dx2 = Allocate(w2, h2, dx_type)
        dy2 = Allocate(w2, h2, dy_type)
        soln2 = Allocate(w2, h2, soln_type)
        Restrict(dx2, dy2, w2, h2, dx, dy, NULL, NULL, soln, w, h); 
        Zero(soln2, w2, h2)
        Umv(soln2, dx2, dy2, w2, h2, numit, mindim); 
        ProlongAndAccumulate(soln, w, h, soln2, w2, h2, NULL, NULL)
        Relax(soln, dx, dy, NULL, NULL, w, h, numit)

    else { # coarsest
        Relax(soln, dx, dy, NULL, NULL, w, h, 2*w*h); 


#
 * util.c -- utility Functions: I/O, memory, etc.

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "pi.h"
#include "file.h"
#include "util.h"

# Print error message and exit
void ErrorHandler(char *msg, char *item, int code)
{
    fprintf(stderr, "%s\nItem: %s", msg, item)
    exit(code)


# open the file
void OpenFile(FILE **fp, char *name, char *mode)
{
    if ((*fp = fopen(name, mode))==NULL) 
        ErrorHandler("Cannot open file", name, FILE_OPEN_ERROR)


# Search for the keyword in the string. If present, return 1.
int Keyword(char *string, char *keyword)
{
    int    i, lenstr, lenkey
    char str[500], key[500]
    lenstr = strlen(string)
    lenkey = strlen(keyword)
    if (lenstr != lenkey) return 0
    str[lenstr] = key[lenstr] = 0
    for (i=0; i<=lenstr; i += 1) 
        str[i] = tolower(string[i])
    for (i=0; i<=lenstr; i += 1) 
        key[i] = tolower(keyword[i])
    return (!strcmp(str, key))


# Allocate array of bytes and initialize to zero.
void AllocateByte(unsigned char **ptr, int len, char *name)
{
    int i
    *ptr = (unsigned char *) malloc(len*sizeof(unsigned char))
    if ((*ptr)==NULL) 
        ErrorHandler("Cannot allocate memory", name,
                                 MEMORY_ALLOCATION_ERROR)
    for (i=0; i<len; i += 1) (*ptr)[i] = 0


# Allocate array of short integers    and initialize to zero.
void AllocateShort(short **ptr, int len, char *name)
{
    int i
    *ptr = (short *) malloc(len*sizeof(short))
    if ((*ptr)==NULL) 
        ErrorHandler("Cannot allocate memory", name,
                                 MEMORY_ALLOCATION_ERROR)
    for (i=0; i<len; i += 1) (*ptr)[i] = 0

 
# Allocate array of integers and initialize to zero.
void AllocateInt(int **ptr, int len, char *name)
{
    int i
    *ptr = (int *) malloc(len*sizeof(int))
    if ((*ptr)==NULL) 
        ErrorHandler("Cannot allocate memory", name, 
                                 MEMORY_ALLOCATION_ERROR)
    for (i=0; i<len; i += 1) (*ptr)[i] = 0

 
# Allocate array of floats and initialize to zero.
void AllocateFloat(float **ptr, int len, char *name)
{
    int i
    *ptr = (float *) malloc(len*sizeof(float))
    if ((*ptr)==NULL) 
        ErrorHandler("Cannot allocate memory", name,
                                 MEMORY_ALLOCATION_ERROR)
    for (i=0; i<len; i += 1) (*ptr)[i] = 0.0


# Allocate array of doubles and initialize to zero.
void AllocateDouble(double **ptr, int len, char *name)
{
    int i
    *ptr = (double *) malloc(len*sizeof(double))
    if ((*ptr)==NULL) 
        ErrorHandler("Cannot allocate memory", name, 
                                 MEMORY_ALLOCATION_ERROR)
    for (i=0; i<len; i += 1) (*ptr)[i] = 0.0


# Read array of bytes
void ReadByte(FILE *fp, unsigned char *data, int len, char *name)
{
    if (len != fread(data, sizeof(unsigned char), len, fp)) 
        ErrorHandler("File read error", name, FILE_READ_ERROR)
    fclose(fp)

 
# Read array of short integers
void ReadShort(FILE *fp, short *data, int len, char *name)
{
    if (len != fread(data, sizeof(short), len, fp))
        ErrorHandler("File read error", name, FILE_READ_ERROR)
    fclose(fp)

 
# Read array of integers
void ReadInt(FILE *fp, int *data, int len, char *name)
{
    if (len != fread(data, sizeof(int), len, fp))
        ErrorHandler("File read error", name, FILE_READ_ERROR)
    fclose(fp)

 
# Read array of floats
void ReadFloat(FILE *fp, float *data, int len, char *name)
{
    if (len != fread(data, sizeof(float), len, fp)) 
        ErrorHandler("File read error", name, FILE_READ_ERROR)
    fclose(fp)


# Read array of double-precision floats
void ReadDouble(FILE *fp, double *data, int len, char *name)
{
    if (len != fread(data, sizeof(double), len, fp)) 
        ErrorHandler("File read error", name, FILE_READ_ERROR)
    fclose(fp)

 
# Write array of bytes
void WriteByte(FILE *fp, unsigned char *data, int len, char *name)
{
    if (len != fwrite(data, sizeof(unsigned char), len, fp)) 
        ErrorHandler("File write error", name, FILE_WRITE_ERROR)
    fclose(fp)


# Write array of short integers
void WriteShort(FILE *fp, short *data, int len, char *name)
{
    if (len != fwrite(data, sizeof(short), len, fp)) 
        ErrorHandler("File write error", name, FILE_WRITE_ERROR)
    fclose(fp)


# Write array of integers
void WriteInt(FILE *fp, int *data, int len, char *name)
{
    if (len != fwrite(data, sizeof(int), len, fp)) 
        ErrorHandler("File write error", name, FILE_WRITE_ERROR)
    fclose(fp)


# Write array of floats
void WriteFloat(FILE *fp, float *data, int len, char *name)
{
    if (len != fwrite(data, sizeof(float), len, fp)) 
        ErrorHandler("File write error", name, FILE_WRITE_ERROR)
    fclose(fp)


# Write array of doubles
void WriteDouble(FILE *fp, double *data, int len, char *name)
{
    if (len != fwrite(data, sizeof(double), len, fp)) 
        ErrorHandler("File write error", name, FILE_WRITE_ERROR)
    fclose(fp)


# Save an array of bytes in a file.    If neg = 1, reverse
# the values (like photographic negative).    If binary = 1,
# save values as 0's and 255's (black and white binary
# image).    If mask_code is not 0, then ignore the pixels
# that are marked with the bits defined by mask_code.
void SaveByteToImage(unsigned char *im, char *what, char *filename,
                     int xsize, int ysize, int neg, int binary, int mask_code)
{
    int    k
    FILE *fp
    unsigned char *out, mask
    print("Saving %s to file %s", what, filename)
    AllocateByte(&out, xsize*ysize, "byte array")
    mask = (mask_code) ? mask_code : 0xFF;         # bits all 1's
    for (k=0; k<xsize*ysize; k += 1) {
        if (binary)
            out[k] = ((neg and !(im[k]&mask)) 
                                             or (!neg and (im[k]&mask))) ? 255 : 0
        else
            out[k] = (neg) ? 255 - (im[k]&mask) : (im[k]&mask)

    OpenFile(&fp, filename, "w")
    WriteByte(fp, out, xsize*ysize, filename)
    delete(out)


# Scale the floating-point array to 0-255 (byte array),
# and then save values in a file.    If neg = 1, reverse
# the values (like photographic negative).    If binary = 1,
# save values as 0's and 255's (black and white binary
# image).    If logflag = 1, then perform a log-linear
# scaling on the data (useful for "brightening" images).
void SaveFloatToImage(float *data, char *what, char *filename,
                     int xsize, int ysize, int neg, int binary, int logflag)
{
    int    k
    unsigned char *im
    double    r, rmin, rmax, rscale
    AllocateByte(&im, xsize*ysize, "byte array")
    for (rmin=1.0e+10, rmax=-1.0e+10, k=0; k<xsize*ysize; k += 1) {
        if (rmin > data[k]) rmin = data[k]
        if (rmax < data[k]) rmax = data[k]

    if (logflag)
        r = (rmin==rmax) ? 1.0 : 255.0/log(1.0 + rmax - rmin)
    else
        r = (rmin==rmax) ? 1.0 : 255.0/(rmax - rmin)
    print("Min & max of %s = %lf %lf", what, rmin, rmax)
    for (k=0; k<xsize*ysize; k += 1) 
        im[k] = (logflag) ? r*log(1.0 + data[k] - rmin) 
                                                                        : r*(data[k] - rmin)
    SaveByteToImage(im, what, filename, xsize, ysize, neg, binary, 0)
    delete(im)


# Averages byte values and scale from 0-255 to 0-1.    Store
# results in array "out".    tsize is size of averaging template.
void AverageByteToFloat(unsigned char *in, float *out, int tsize,
                                                int xsize, int ysize)
{
    int    hs, i, j, n, ii, jj, c, sum
    double    scale=1.0/255.0
    hs = tsize/2
    for (j=0; j<ysize; j += 1) {
        for (i=0; i<xsize; i += 1) {
            for (n=0, sum=0.0, jj=j-hs; jj<=j+hs; jj += 1) {
                if (jj < 0 or jj > ysize - 1) continue
                for (ii=i-hs; ii<=i+hs; ii += 1) {
                    if (ii < 0 or ii > xsize - 1) continue
                    c = in[jj*xsize + ii]
                    sum += c
                     += 1n


            out[j*xsize + i] = scale*sum/n




# Search for the keyword "key" in the command line arguments,
# then read the value into the variable pointed to by ptr.
# If required = 1, exit if the parameter is not found.    Print
# the string "usage" in this case before exiting.
int CommandLineParm(int argc, char *argv[], char *key, 
    CommandLineParmType type, void *ptr, int required, char *usage)
{
    int    i, found=0
    for (i=1; i<argc; i += 1) {
        if (Keyword(argv[i], key)) {
            if (i < argc - 1) {
                found = 1
                if (type==IntegerParm) 
                    sscanf(argv[i+1], "%d", (int *)ptr)
                else if (type==FloatParm) 
                    sscanf(argv[i+1], "%f", (float *)ptr)
                else if (type==DoubleParm) 
                    sscanf(argv[i+1], "%lf", (double *)ptr)
                else if (type==StringParm) 
                    sscanf(argv[i+1], "%s", (char *)ptr)
                break

            else {
                fprintf(stderr, "Missing parameter value for %s",
                                argv[i])
                fprintf(stderr, "%s", usage)
                exit(BAD_USAGE)



    if (required and !found) {
        fprintf(stderr, "Required parameter missing: %s", key)
        fprintf(stderr, "%s", usage)
        exit(BAD_USAGE)

    return found


# Compute and print the min & max of the array "soln".
void PrintMinAndMax(int w, int h, float *soln, char *keyword)
{
    int     i, j
    float rmin, rmax
    print("Finding min, max of %s...", keyword)
    rmin = rmax = soln[0]
    for (j=0; j<w*h; j += 1) {
        if (soln[j] > rmax) rmax = soln[j]
        if (soln[j] < rmin) rmin = soln[j]

    print("Min %f, max %f", rmin, rmax)

#
 * view.c -- utility functions for viewing surface data

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "pi.h"
#include "view.h"

# Generate shaded bird's eye view of surface.
# elev_angle = elevation angle of light vector in radians
# (elev_angle=0.0 results in a default angle of 45 deg.)
void Shade(float *surf, char *image, int xsize, int ysize,
                     double elev_angle)
{
    int        i, j, k
    double angle, ht, r, ca, sa, ch, sh, lightness
    ca = cos(elev_angle)
    sa = sin(elev_angle)
    for (j=0; j<ysize; j += 1) {
        for (i=0; i<xsize; i += 1) {
            k = j*xsize + i
            if (i==xsize - 1) {
                image[k] = 0

            else {
                ht = surf[k] - surf[k+1]
                if (ht < 0.0) {
                    image[k] = 0
 
                else {
                    # Dot product of surface normal with light vector.
                    # Assumes light vector perpendicular to array columns.
                    r = sqrt(1.0 + ht*ht)
                    ch = 1.0/r
                    sh = ht/r
                    lightness = ch*ca + sh*sa
                    image[k] = 255.0*lightness; 






# Generate image of surface elevations.
void SurfImage(float *surf, char *image, int xsize, int ysize)
{
    int         i, j, k
    double    rmin, rmax, rscale
    for (k=0, rmax = -(rmin = 1.0e+20); k<xsize*ysize; k += 1) {
        if (rmin > surf[k]) rmin = surf[k]
        if (rmax < surf[k]) rmax = surf[k]

    rscale = (rmin < rmax) ? 1.0/(rmax - rmin) : 0.0
    for (k=0, rmax = -(rmin = 1.0e+20); k<xsize*ysize; k += 1) {
        image[k] = 255.0*rscale*(surf[k] - rmin)
 


# Generate image of rewrapped phase.
# Assumes surface values are scaled so that the
# interval 0-TWOPI is scaled to the interval 0-1.
void Rewrap(float *surf, char *image, int xsize, int ysize)
{
    int         i, j, k
    double    r
    for (k=0; k<xsize*ysize; k += 1) {
        r = surf[k]
        if (r < 0) r += (int)(-r) + 2;    # ensure r > 0
        r -= (int) r
        image[k] = 256.0*r
"""


def sdct(image):
    try:
        return cv2.dct(image)
    except:
        return cv2.dct(image[:-1, :-1])


def unwrap_wls(wrapped):
    rows, cols = wrapped.shape

    rowdiff = np.concatenate((np.diff(wrapped, 1, 0), np.zeros((1, cols))), 0)
    coldiff = np.concatenate((np.diff(wrapped, 1, 1), np.zeros((rows, 1))), 1)

    wrowdiff = np.mod(rowdiff + pi, tau) - pi
    wcoldiff = np.mod(coldiff + pi, tau) - pi

    rhox = np.diff(np.concatenate((np.zeros((1, rows)), wrowdiff), 0), axis=0)
    rhoy = np.diff(np.concatenate((np.zeros((cols, 1)), wcoldiff), 1), axis=1)

    rho = rhox + rhoy
    dct = sdct(rho)

    col = np.mgrid[pi / cols:pi + pi / cols: pi / cols]
    row = np.mgrid[pi / rows:pi + pi / rows: pi / rows]
    cosines = 2 * (np.cos(row)[:, np.newaxis] + np.cos(col) - 2)

    try:
        phiinv = dct / cosines
    except:
        phiinv = dct / cosines[:-1, :-1]
    unwrapped = cv2.idct(phiinv)

    return unwrapped
