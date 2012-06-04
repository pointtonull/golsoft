#!/usr/bin/env python
#-*- coding: UTF-8 -*-

"""
Functions for branch cutting
============================
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

    if (a == c and b == d):
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
                        best_dist2 = dist2;
                        besta = i;
                        bestb = j;
        if found:
            ra = besta
            rb = bestb
            break
    return best_dist2;

"""
Functions for finding dipoles and eliminating them
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
Extract phase from input data
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
                x *= r;
                y *= r;
                angle = atan2((double)y, (double)x);
                if (angle  0) angle += TWOPI;
                if (angle >= TWOPI) angle -= TWOPI;
                angle *= one_over_twopi;
                phase[j*xsize + i] = angle;
            }
        }
    }
    else if (in_format==2) {     /* quantized phase */
        scale = 1.0/256.0;
        for (j=0; j<xsize*ysize; j+= 1) {
            phase[j] = quantized_phase[j]*scale;
        }
    }
    else {    /* 4-byte float phase */
        scale = one_over_twopi;
        for (j=0; j<xsize*ysize; j+= 1) {
            /* re-scale phase to interval (0,1) */
            phase[j] = float_phase[j]*scale;
        }
    }

    /*
     * gold.c -- generate branch cuts by Goldstein's algorithm
     */
    #include <stdio.h>
    #include <math.h>
    #include "util.h"
    #include "brcut.h"
    #include "list.h"

    /* Goldstein's phase-unwrapping algorithm.  The bitflags store */
    /* the masked pixels (to be ignored) and the residues and      */
    /* accumulates other info such as the branch cut pixels.       */
    void GoldsteinBranchCuts(unsigned char *bitflags, int MaxCutLen,
                   int NumRes, int xsize, int ysize, int branchcut_code)
    {
      int            i, j, k, ii, jj, kk, m, n, ri, rj;
      int            charge, boxctr_i, boxctr_j, boxsize, bs2;
      int            dist, min_dist, rim_i, rim_j, near_i, near_j;
      int            ka, num_active, max_active, *active_list;
      int            bench;
      int            draw_cut_line;
      double         r;
      bench = ysize/100;
      if (bench  1) bench = 1;
      if (MaxCutLen  2) MaxCutLen = 2;
      max_active = NumRes + 10;
      AllocateInt(&active_list, max_active + 1, "book keeping data");
      /* branch cuts */
      printf("Computing branch cuts\n");
      for (j=0; j<ysize; j+= 1) {
        if (j%bench==0) {
          printf("%d ", j/bench);
          fflush(stdout);
        }
        for (i=0; i<xsize; i+= 1) {
          k = j*xsize + i;
          if ((bitflags[k] & (POS_RES | NEG_RES))
                   && !(bitflags[k] & VISITED)) {
            bitflags[k] |= VISITED;  /* turn on visited flag */
            bitflags[k] |= ACTIVE;   /* turn on active flag */
            charge = (bitflags[k] & POS_RES) ? 1 : -1;
            num_active = 0;
            active_list[num_active+= 1] = k;
            if (num_active > max_active) num_active = max_active;
            for (boxsize = 3; boxsize=2*MaxCutLen; boxsize += 2) {
              bs2 = boxsize/2;
              for (ka=0; ka<num_active; ka+= 1) {
                boxctr_i = active_list[ka]%xsize;
                boxctr_j = active_list[ka]/xsize;
                for (jj=boxctr_j - bs2; jj=boxctr_j + bs2; jj+= 1) {
                  for (ii=boxctr_i - bs2; ii=boxctr_i + bs2; ii+= 1) {
                    kk = jj*xsize + ii;
                    if (ii<0 || ii>=xsize || jj<0 || jj>=ysize) {
                      continue;
                    }
                    else {
                      if (ii==0 || ii==xsize-1 || jj==0 || jj==ysize-1
                            || (bitflags[kk] & BORDER)) {
                        charge = 0;
                        DistToBorder(bitflags, BORDER, boxctr_i,
                                 boxctr_j, &ri, &rj, xsize, ysize);
                        PlaceCut(bitflags, ri, rj, boxctr_i, boxctr_j,
                                 xsize, ysize, branchcut_code);
                      }
                      else if ((bitflags[kk] & (POS_RES | NEG_RES))
                                       && !(bitflags[kk] & ACTIVE)) {
                        if (!(bitflags[kk] & VISITED)) {
                          charge += (bitflags[kk] & POS_RES) ? 1 : -1;
                          bitflags[kk] |= VISITED;   /* set flag */
                        }
                        active_list[num_active+= 1] = kk;
                        if (num_active > max_active)
                               num_active = max_active;
                        bitflags[kk] |= ACTIVE;  /* set active flag */
                        PlaceCut(bitflags, ii, jj, boxctr_i, boxctr_j,
                                 xsize, ysize, branchcut_code);
                      }
                      if (charge==0)
                        goto continue_scan;
                    }  /* else */
                  }   /* for (ii ... */
                }   /* for (jj ... */
              }  /* for (ka ... */
            }   /* for (boxsize ... */

            if (charge != 0) {   /* connect branch cuts to rim */
              min_dist = xsize + ysize;  /* large value */
              for (ka=0; ka<num_active; ka+= 1) {
                ii = active_list[ka]%xsize;
                jj = active_list[ka]/xsize;
                if ((dist = DistToBorder(bitflags, BORDER,
                            ii, jj, &ri, &rj, xsize, ysize))<min_dist) {
                  min_dist = dist;
                  near_i = ii;
                  near_j = jj;
                  rim_i = ri;
                  rim_j = rj;
                }
              }
              PlaceCut(bitflags, near_i, near_j, rim_i, rim_j,
                       xsize, ysize, branchcut_code);
            }
            continue_scan :
            /* mark all active pixels inactive */
            for (ka=0; ka<num_active; ka+= 1)
              bitflags[active_list[ka]] &= ~ACTIVE;  /* turn flag off */
          }  /* if (bitflags ... */
        }  /* for (i ... */
      }  /* for (j ... */
      printf("\n");
      free(active_list);
      return;
    }
    /*
     *   grad.c -- function for computing phase derivative
     *             (wrapped phase difference)
     */
    #include <stdio.h>
    #include <math.h>
    #include "grad.h"
    /* return wrapped phase difference */
    float Gradient(float p1, float p2)
    {
      float  r;
      r = p1 - p2;
      if (r > 0.5) r -= 1.0;
      if (r  -0.5) r += 1.0;
      return r;
    }
    /*
     * maskfat.c -- fatten mask
     */
    #include <stdlib.h>
    #include <stdio.h>
    #include <string.h>
    #include <math.h>
    #include "util.h"

    /* Fattens mask and border by "thickness" pixels */
    void FattenMask(unsigned char *array, int mask_code, int thickness,
                    int xsize, int ysize)
    {
      int i, j, k, ii, jj, kk, t2=thickness;
      unsigned char *temp;
      if (t2  1) return;
      AllocateByte(&temp, xsize*ysize, "temp array");
      for (j=0; j<ysize; j+= 1) {
        for (i=0; i<xsize; i+= 1) {
          k = j*xsize + i;
          temp[k] = 0;
          for (jj=j-t2; jj=j+t2 && !temp[k]; jj+= 1) {
            for (ii=i-t2; ii=i+t2 && !temp[k]; ii+= 1) {
              kk = jj*xsize + ii;
              if (ii>=0 && ii<xsize && jj>=0 && jj<ysize
                                         && (array[kk] & mask_code))
                temp[k] = 1;
            }
          }
        }
      }
      for (k=0; k<xsize*ysize; k+= 1) {
        if (temp[k]) array[k] |= mask_code;
      }
      free(temp);
    }

    /* Fattens zero values in quality map by "thickness" pixels */
    void FattenQual(float *qual_map, int thickness,
                    int xsize, int ysize)
    {
      int i, j, k, ii, jj, kk, t2=thickness;
      unsigned char *temp;
      if (t2  1) return;
      AllocateByte(&temp, xsize*ysize, "temp array");
      for (j=0; j<ysize; j+= 1) {
        for (i=0; i<xsize; i+= 1) {
          k = j*xsize + i;
          temp[k] = 0;
          for (jj=j-t2; jj=j+t2 && !temp[k]; jj+= 1) {
            for (ii=i-t2; ii=i+t2 && !temp[k]; ii+= 1) {
              kk = jj*xsize + ii;
              if (ii>=0 && ii<xsize && jj>=0 && jj<ysize
                          && qual_map[kk]==0.0)
                temp[k] = 1;
            }
          }
        }
      }
      for (k=0; k<xsize*ysize; k+= 1) {
        if (temp[k]) qual_map[k] = 0.0;
      }
      free(temp);
    }
/*
 * maingold.c - phase unwrapping by means of residues & branch cuts
 *
 * Source code files required:
 *        brcut.c      dipole.c     extract.c        gold.c
 *         grad.c        list.c    maingold.c     maskfat.c
 *         path.c    residues.c        util.c
 */
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
  int            i, j, k, m, n, a, b;
  FILE           *ifp, *ofp, *mfp;
  float          *phase;     /* array */ 
  float          *soln;      /* array */
  unsigned char  *bitflags;
  char           buffer[200], string[200];
  char           infile[200], outfile[200];
  char           maskfile[200], tempstr[200], format[200];
  int            in_format, debug_flag, dipole_flag;
  int            xsize, ysize;   /* dimensions of arrays */ 
  int            NumRes, MaxCutLen;
  char           use[] =    /* define usage statement */
    "Usage: prog-name -input file -format fkey -output file\n"
    "  -xsize x -ysize y [ -mask file -cutlen c -debug yes/no\n"
    "  -dipole yes/no ]\n"
    "where 'fkey' is a keyword designating the input file type\n"
    "(key = complex8, complex4, float or byte), 'x' and 'y' are\n"
    "the dimensions of the file, mask is an optional byte-file\n"
    "of masks for masking out undefined phase values and cutlen\n"
    "is the max branchcut length allowed.  All files are simple\n"
    "raster files, and the output file consists of floating\n"
    "point numbers that define the heights of the unwrapped\n"
    "surface.  If the 'debug' parm is 'yes', then intermediate\n"
    "byte-files are saved (residues, branch cuts, etc.)  If\n"
    "the 'dipole' parm is 'yes', then dipole-residues are\n"
    "eliminated before unwrapping.\n";

  printf("Phase Unwrapping By Means of Goldstein's Algorithm\n");

  /* GET COMMAND LINE PARAMETERS AND CHECK */
  CommandLineParm(argc,argv, "-input", StringParm, infile, 1, use);
  CommandLineParm(argc,argv, "-format", StringParm, format, 1, use);
  CommandLineParm(argc,argv, "-output", StringParm, outfile, 1,use);
  CommandLineParm(argc, argv, "-xsize", IntegerParm, &xsize, 1,use);
  CommandLineParm(argc, argv, "-ysize", IntegerParm, &ysize, 1,use);
  if (!CommandLineParm(argc, argv, "-mask", StringParm, maskfile,
        0, use)) strcpy(maskfile, "none");
  if (!CommandLineParm(argc, argv, "-cutlen", IntegerParm,
     &MaxCutLen, 0, use)) MaxCutLen = 0; /* default defined below */
  if (!CommandLineParm(argc, argv, "-debug", StringParm, tempstr,
        0, use)) debug_flag = 0;
  else debug_flag = Keyword(tempstr, "yes");
  if (!CommandLineParm(argc, argv, "-dipole", StringParm, tempstr,
        0, use)) dipole_flag=0;
  else dipole_flag = Keyword(tempstr, "yes");

  if (Keyword(format, "complex8"))  in_format = 0;
  else if (Keyword(format, "complex4"))  in_format = 1;
  else if (Keyword(format, "byte"))  in_format = 2;
  else if (Keyword(format, "float"))  in_format = 3;
  else {
    fprintf(stderr, "Unrecognized format: %s\n", format);
    exit(BAD_PARAMETER);
  }

  printf("Input file =  %s\n", infile);
  printf("Input file type = %s\n", format);
  printf("Output file =  %s\n", outfile);
  printf("File dimensions = %dx%d (cols x rows).\n", xsize, ysize);

  if (Keyword(maskfile, "none")) printf("No mask file.\n");
  else printf("Mask file = %s\n", maskfile);

  if (dipole_flag) printf("Dipole-residues will be eliminated.\n");
  if (debug_flag) 
    printf("Intermediate files will be saved (i.e. debug on).\n");

  /*  OPEN FILES, ALLOCATE MEMORY      */
  OpenFile(&ifp, infile, "rb");
  OpenFile(&ofp, outfile, "wb");
  mfp = NULL;
  if (!Keyword(maskfile, "none"))
    OpenFile(&mfp, maskfile, "rb");

  AllocateFloat(&phase, xsize*ysize, "phase data");
  AllocateFloat(&soln, xsize*ysize, "unwrapped data");
  AllocateByte(&bitflags, xsize*ysize, "bitflag data");

  /*  READ AND PROCESS DATA  */
  printf("Reading phase data...\n");
  GetPhase(in_format, ifp, infile, phase, xsize, ysize);

  /* mask data */
  printf("Processing mask data...\n");
  if (mfp) {
    ReadByte(mfp, bitflags, xsize*ysize, maskfile);
  }
  else {
    for (k=0; k<xsize*ysize; k+= 1)
      bitflags[k] = 255;
  }
  for (k=0; k<xsize*ysize; k+= 1)
    bitflags[k] = (!bitflags[k]) ? BORDER : 0;

  FattenMask(bitflags, BORDER, 1, xsize, ysize);

  /*  LOCATE AND PROCESS RESIDUES  */
  /* compute residues and store in bitflags array */
  NumRes = Residues(phase, bitflags, POS_RES, NEG_RES,
                    BORDER, xsize, ysize);
  printf("%d Residues\n", NumRes);

  if (debug_flag) {
    char filename[300];
    sprintf(filename, "%s.res", outfile);
    SaveByteToImage(bitflags, "residues", filename, xsize, ysize,
                    1, 1, 0);
  }

  /*  GENERATE BRANCH CUTS  */
  if (dipole_flag) {  /* elimate dipole-residues first */
    Dipole(bitflags, xsize, ysize, BRANCH_CUT);
    i = Residues(phase, bitflags, 0, 0, BORDER | BRANCH_CUT,
                 xsize, ysize);
    printf("%d Residues are left\n", i);
  }

  if (MaxCutLen==0) MaxCutLen = (xsize + ysize)/2;
  GoldsteinBranchCuts(bitflags, MaxCutLen, NumRes, xsize, ysize,
                      BRANCH_CUT);
  if (debug_flag) {
    char filename[300];
    sprintf(filename, "%s.brc", outfile);
    SaveByteToImage(bitflags, "branch cuts", filename, 
                    xsize, ysize, 1, 1, BRANCH_CUT | BORDER);
  }

  /*  UNWRAP AROUND CUTS */
  for (k=0, i=0; i<xsize*ysize; i+= 1) 
    if (bitflags[i] & BRANCH_CUT) k+= 1;
  printf("%d BRANCH CUT PIXELS\n", k);
  printf("Unwrapping around branch cuts\n");
  k = UnwrapAroundCuts(phase, bitflags, soln, xsize, ysize, AVOID,
                       0, NULL);
  if (k > 1) printf("%d disconnected pieces.\n", k);
  else printf("%d piece\n", k);
  printf("\nFinished\n");
  PrintMinAndMax(xsize, ysize, soln, "solution");

  /*  SAVE RESULT  */
  for (k=0; k<xsize*ysize; k+= 1)
    soln[k] *= TWOPI; /* scale output */
  printf("Saving unwrapped surface to file '%s'\n", outfile);
  WriteFloat(ofp, soln, xsize*ysize, outfile);
  free(soln);
  free(phase);
  free(bitflags);
}
    /*
     *   path.c -- function for unwrapping around cuts
     */
    #include <stdio.h>
    #include <math.h>
    #include "list.h"
    #include "grad.h"
    //#include "util.h"
    #include "path.h"

    /* Unwrap the phase data (by Itoh's method) without crossing
     * any branch cuts.  Return number of disconnected pieces.
     */
    int UnwrapAroundCuts(float *phase, unsigned char *bitflags,
                  float *soln, int xsize, int ysize, int cut_code,
                  int debug_flag, char *infile)
    {
      int  i, j, k, a, b, c, n, num_pieces=0;
      float  value;
      float  min_qual, small_val = -1.0E+10;
      int    num_index, max_list_size, bench, benchout;
      int    *index_list;
      int    unwrapped_code=UNWRAPPED, postponed_code = POSTPONED;
      int    avoid_code;

      bench = xsize*ysize/100;
      if (bench  1) bench = 1;
      benchout = 10*bench;
      min_qual = small_val;
      max_list_size = xsize*ysize; /* this size may be reduced */
      AllocateInt(&index_list, max_list_size + 1,
                  "bookkeeping list (index)");

      avoid_code = cut_code | unwrapped_code | BORDER;

      /* find starting point */
      n = 0;
      num_index = 0;
      for (j=0; j<ysize; j+= 1) {
        for (i=0; i<xsize; i+= 1) {
          k = j*xsize + i;
          if (!(bitflags[k] & avoid_code)) {
            bitflags[k] |= unwrapped_code;
            if (bitflags[k] & postponed_code)
              /* soln[k] already stores the unwrapped value */
              value = soln[k];
            else {
              += 1num_pieces;
              value = soln[k] = phase[k];
            }
            UpdateList(NULL, i, j, value, phase, soln, bitflags, xsize,
                       ysize, index_list, &num_index, avoid_code,
                       unwrapped_code, postponed_code, max_list_size,
                       0, &min_qual);
            while (num_index > 0) {
              if (n%bench==0) {
                printf("%d ", n/bench);
                fflush(stdout);
                if (0 && debug_flag && n%benchout==0 && n>0) {
                  char filename[300];
                  static char file_num = 0;
                  sprintf(filename, "%s.fill-%d", infile, += 1file_num);
                  SaveByteToImage(bitflags, "fill path", filename,
                                  xsize, ysize, 1, 1, avoid_code);
                }
              }
              += 1n;
              if (!GetNextOneToUnwrap(&a, &b, index_list,
                                      &num_index, xsize, ysize))
                break;   /* no more to unwrap */
              c = b*xsize + a;
              bitflags[c] |= unwrapped_code;
              value = soln[c];
              UpdateList(NULL, a, b, value, phase, soln, bitflags,
                         xsize, ysize, index_list, &num_index,
                         avoid_code, unwrapped_code, postponed_code,
                         max_list_size, 0, &min_qual);
            }
          }
        }
      }
      printf("\n");
      free(index_list);
      /* unwrap branch cut pixels */
      for (j=1; j<ysize; j+= 1) {
        for (i=1; i<xsize; i+= 1) {
          k = j*xsize + i;
          if (bitflags[k] & cut_code) {
            if (!(bitflags[k-1] & cut_code))
              soln[k] = soln[k-1] + Gradient(phase[k], phase[k-1]);
            else if (!(bitflags[k-xsize] & cut_code))
              soln[k]
                 = soln[k-xsize] + Gradient(phase[k], phase[k-xsize]);
          }
        }
      }
      return num_pieces;
    }
/*
 *   residues.c -- function for extracting residues from phase
 */
#include <stdio.h>
#include <math.h>
#include "residues.h"
#include "grad.h"
/* Detect residues in phase data and mark them as positive or  */
/* negative residues in the bitflags array.  Ignore the pixels */
/* marked with avoid_code in the bitflags araay.               */
int Residues(float *phase, unsigned char *bitflags, int posres_code,
             int negres_code, int avoid_code, int xsize, int ysize)
{
  int  i, j, k, NumRes=0;
  double  r;
  for (j=0; j<ysize - 1; j+= 1) {
    for (i=0; i<xsize - 1; i+= 1) {
      k = j*xsize + i;
      if (bitflags && ((bitflags[k] & avoid_code)
           || (bitflags[k+1] & avoid_code)
             || (bitflags[k+1+xsize] & avoid_code)
                 || (bitflags[k+xsize] & avoid_code))) {
        continue; /* masked region: don't unwrap */
      }
      r = Gradient(phase[k+1], phase[k])
           + Gradient(phase[k+1+xsize], phase[k+1])
              + Gradient(phase[k+xsize], phase[k+1+xsize])
                 + Gradient(phase[k], phase[k+xsize]);
      if (bitflags) {
        if (r > 0.01) bitflags[k] |= posres_code;
        else if (r  -0.01) bitflags[k] |= negres_code;
      }
      if (r*r > 0.01) += 1NumRes;
    }
  }
  return NumRes;
}
    /*
     * util.c -- utility Functions: I/O, memory, etc.
     */
    #include <stdlib.h>
    #include <stdio.h>
    #include <string.h>
    #include <math.h>
    #include "pi.h"
    #include "file.h"
    #include "util.h"

    /* Print error message and exit */
    void ErrorHandler(char *msg, char *item, int code)
    {
      fprintf(stderr, "%s\nItem: %s\n", msg, item);
      exit(code);
    }

    /* open the file */
    void OpenFile(FILE **fp, char *name, char *mode)
    {
      if ((*fp = fopen(name, mode))==NULL)
        ErrorHandler("Cannot open file", name, FILE_OPEN_ERROR);
    }

    /* Search for the keyword in the string. If present, return 1. */
    int Keyword(char *string, char *keyword)
    {
      int  i, lenstr, lenkey;
      char str[500], key[500];
      lenstr = strlen(string);
      lenkey = strlen(keyword);
      if (lenstr != lenkey) return 0;
      str[lenstr] = key[lenstr] = 0;
      for (i=0; i=lenstr; i+= 1)
        str[i] = tolower(string[i]);
      for (i=0; i=lenstr; i+= 1)
        key[i] = tolower(keyword[i]);
      return (!strcmp(str, key));
    }

    /* Allocate array of bytes and initialize to zero. */
    void AllocateByte(unsigned char **ptr, int len, char *name)
    {
      int i;
      *ptr = (unsigned char *) malloc(len*sizeof(unsigned char));
      if ((*ptr)==NULL)
        ErrorHandler("Cannot allocate memory", name,
                     MEMORY_ALLOCATION_ERROR);
      for (i=0; i<len; i+= 1) (*ptr)[i] = 0;
    }

    /* Allocate array of short integers  and initialize to zero. */
    void AllocateShort(short **ptr, int len, char *name)
    {
      int i;
      *ptr = (short *) malloc(len*sizeof(short));
      if ((*ptr)==NULL)
        ErrorHandler("Cannot allocate memory", name,
                     MEMORY_ALLOCATION_ERROR);
      for (i=0; i<len; i+= 1) (*ptr)[i] = 0;
    }

    /* Allocate array of integers and initialize to zero. */
    void AllocateInt(int **ptr, int len, char *name)
    {
      int i;
      *ptr = (int *) malloc(len*sizeof(int));
      if ((*ptr)==NULL)
        ErrorHandler("Cannot allocate memory", name,
                     MEMORY_ALLOCATION_ERROR);
      for (i=0; i<len; i+= 1) (*ptr)[i] = 0;
    }

    /* Allocate array of floats and initialize to zero. */
    void AllocateFloat(float **ptr, int len, char *name)
    {
      int i;
      *ptr = (float *) malloc(len*sizeof(float));
      if ((*ptr)==NULL)
        ErrorHandler("Cannot allocate memory", name,
                     MEMORY_ALLOCATION_ERROR);
      for (i=0; i<len; i+= 1) (*ptr)[i] = 0.0;
    }

    /* Allocate array of doubles and initialize to zero. */
    void AllocateDouble(double **ptr, int len, char *name)
    {
      int i;
      *ptr = (double *) malloc(len*sizeof(double));
      if ((*ptr)==NULL)
        ErrorHandler("Cannot allocate memory", name,
                     MEMORY_ALLOCATION_ERROR);
      for (i=0; i<len; i+= 1) (*ptr)[i] = 0.0;
    }

    /* Read array of bytes */
    void ReadByte(FILE *fp, unsigned char *data, int len, char *name)
    {
     if (len != fread(data, sizeof(unsigned char), len, fp))
            ErrorHandler("File read error", name, FILE_READ_ERROR);
      fclose(fp);
    }

    /* Read array of short integers */
    void ReadShort(FILE *fp, short *data, int len, char *name)
    {
      if (len != fread(data, sizeof(short), len, fp))
        ErrorHandler("File read error", name, FILE_READ_ERROR);
      fclose(fp);
    }

    /* Read array of integers */
    void ReadInt(FILE *fp, int *data, int len, char *name)
    {
      if (len != fread(data, sizeof(int), len, fp))
        ErrorHandler("File read error", name, FILE_READ_ERROR);
      fclose(fp);
    }

    /* Read array of floats */
    void ReadFloat(FILE *fp, float *data, int len, char *name)
    {
      if (len != fread(data, sizeof(float), len, fp))
        ErrorHandler("File read error", name, FILE_READ_ERROR);
      fclose(fp);
    }

    /* Read array of double-precision floats */
    void ReadDouble(FILE *fp, double *data, int len, char *name)
    {
      if (len != fread(data, sizeof(double), len, fp))
        ErrorHandler("File read error", name, FILE_READ_ERROR);
      fclose(fp);
    }

    /* Write array of bytes */
    void WriteByte(FILE *fp, unsigned char *data, int len, char *name)
    {
      if (len != fwrite(data, sizeof(unsigned char), len, fp))
        ErrorHandler("File write error", name, FILE_WRITE_ERROR);
      fclose(fp);
    }

    /* Write array of short integers */
    void WriteShort(FILE *fp, short *data, int len, char *name)
    {
      if (len != fwrite(data, sizeof(short), len, fp))
        ErrorHandler("File write error", name, FILE_WRITE_ERROR);
      fclose(fp);
    }

    /* Write array of integers */
    void WriteInt(FILE *fp, int *data, int len, char *name)
    {
      if (len != fwrite(data, sizeof(int), len, fp))
        ErrorHandler("File write error", name, FILE_WRITE_ERROR);
      fclose(fp);
    }

    /* Write array of floats */    // Modified by GX LIU on Sept. 18, 2001
    void WriteFloat(FILE *fp, float *data, int len, char *name)
    {
      /* int i, num = 0;

      fprintf(fp, "%s\n", "DSAA");
      fprintf(fp, "%d %d\n", 4872, 5296);
      fprintf(fp, "%d %d\n", 0, 4871);
      fprintf(fp, "%d %d\n", 0, 5295);
      fprintf(fp, "%f %f\n", -80.00, 200.000);

      for (i=0; i<len; i+= 1) {
        num += fprintf(fp, "%8.2f ", data[i]);
        if (num==4872) {
            fprintf(fp, "\n"); num=0; }
      } */

      if (len != fwrite(data, sizeof(float), len, fp))
            ErrorHandler("File write error", name, FILE_WRITE_ERROR);
      fclose(fp);

    }

    /* Write array of doubles */
    void WriteDouble(FILE *fp, double *data, int len, char *name)
    {
      if (len != fwrite(data, sizeof(double), len, fp))
        ErrorHandler("File write error", name, FILE_WRITE_ERROR);
      fclose(fp);
    }

    /* Save an array of bytes in a file.  If neg = 1, reverse   */
    /* the values (like photographic negative).  If binary = 1, */
    /* save values as 0's and 255's (black and white binary     */
    /* image).  If mask_code is not 0, then ignore the pixels   */
    /* that are marked with the bits defined by mask_code.      */
    void SaveByteToImage(unsigned char *im, char *what, char *filename,
               int xsize, int ysize, int neg, int binary, int mask_code)
    {
      int  k;
      FILE *fp;
      unsigned char *out, mask;
      printf("Saving %s to file %s\n", what, filename);
      AllocateByte(&out, xsize*ysize, "byte array");
      mask = (mask_code) ? mask_code : 0xFF;     /* bits all 1's */
      for (k=0; k<xsize*ysize; k+= 1) {
        if (binary)
          out[k] = ((neg && !(im[k]&mask))
                           || (!neg && (im[k]&mask))) ? 255 : 0;
        else
          out[k] = (neg) ? 255 - (im[k]&mask) : (im[k]&mask);
      }
      OpenFile(&fp, filename, "w");
      WriteByte(fp, out, xsize*ysize, filename);
      free(out);
    }

    /* Scale the floating-point array to 0-255 (byte array),    */
    /* and then save values in a file.  If neg = 1, reverse     */
    /* the values (like photographic negative).  If binary = 1, */
    /* save values as 0's and 255's (black and white binary     */
    /* image).  If logflag = 1, then perform a log-linear       */
    /* scaling on the data (useful for "brightening" images).   */
    void SaveFloatToImage(float *data, char *what, char *filename,
               int xsize, int ysize, int neg, int binary, int logflag)
    {
      int  k;
      unsigned char *im;
      double  r, rmin, rmax, rscale;
      AllocateByte(&im, xsize*ysize, "byte array");
      for (rmin=1.0e+10, rmax=-1.0e+10, k=0; k<xsize*ysize; k+= 1) {
        if (rmin > data[k]) rmin = data[k];
        if (rmax  data[k]) rmax = data[k];
      }
      if (logflag)
        r = (rmin==rmax) ? 1.0 : 255.0/log(1.0 + rmax - rmin);
      else
        r = (rmin==rmax) ? 1.0 : 255.0/(rmax - rmin);
      printf("Min & max of %s = %lf %lf\n", what, rmin, rmax);
      for (k=0; k<xsize*ysize; k+= 1)
        im[k] = (logflag) ? r*log(1.0 + data[k] - rmin)
                                        : r*(data[k] - rmin);
      SaveByteToImage(im, what, filename, xsize, ysize, neg, binary, 0);
      free(im);
    }

    /* Averages byte values and scale from 0-255 to 0-1.  Store      */
    /* results in array "out".  tsize is size of averaging template. */
    void AverageByteToFloat(unsigned char *in, float *out, int tsize,
                            int xsize, int ysize)
    {
      int  hs, i, j, n, ii, jj, c, sum;
      double  scale=1.0/255.0;
      hs = tsize/2;
      for (j=0; j<ysize; j+= 1) {
        for (i=0; i<xsize; i+= 1) {
          for (n=0, sum=0.0, jj=j-hs; jj=j+hs; jj+= 1) {
            if (jj  0 || jj > ysize - 1) continue;
            for (ii=i-hs; ii=i+hs; ii+= 1) {
              if (ii  0 || ii > xsize - 1) continue;
              c = in[jj*xsize + ii];
              sum += c;
              += 1n;
            }
          }
          out[j*xsize + i] = scale*sum/n;
        }
      }
    }

    /* Search for the keyword "key" in the command line arguments, */
    /* then read the value into the variable pointed to by ptr.    */
    /* If required = 1, exit if the parameter is not found.  Print */
    /* the string "usage" in this case before exiting.             */
    int CommandLineParm(int argc, char *argv[], char *key,
      CommandLineParmType type, void *ptr, int required, char *usage)
    {
      int  i, found=0;
      for (i=1; i<argc; i+= 1) {
        if (Keyword(argv[i], key)) {
          if (i  argc - 1) {
            found = 1;
            if (type==IntegerParm)
              sscanf(argv[i+1], "%d", (int *)ptr);
            else if (type==FloatParm)
              sscanf(argv[i+1], "%f", (float *)ptr);
            else if (type==DoubleParm)
              sscanf(argv[i+1], "%lf", (double *)ptr);
            else if (type==StringParm)
              sscanf(argv[i+1], "%s", (char *)ptr);
            break;
          }
          else {
            fprintf(stderr, "Missing parameter value for %s\n",
                    argv[i]);
            fprintf(stderr, "%s", usage);
            exit(BAD_USAGE);
          }
        }
      }
      if (required && !found) {
        fprintf(stderr, "Required parameter missing: %s\n", key);
        fprintf(stderr, "%s", usage);
        exit(BAD_USAGE);
      }
      return found;
    }

    /* Compute and print the min & max of the array "soln". */
    void PrintMinAndMax(int w, int h, float *soln, char *keyword)
    {
      int   i, j;
      float rmin, rmax;
      printf("Finding min, max of %s...", keyword);
      rmin = rmax = soln[0];
      for (j=0; j<w*h; j+= 1) {
        if (soln[j] > rmax) rmax = soln[j];
        if (soln[j]  rmin) rmin = soln[j];
      }
      printf("Min %f, max %f\n", rmin, rmax);
    }


def main():
    if len(sys.argv) == 1:
        print("Se debe especificar al menos un fichero de entrada.")
    else:
        image = Image.open(sys.argv[1])
        uw_image = unwrap_image(image)

        if len(sys.argv) == 2:
            uw_image.show()
            return True
        elif len(sys.argv) == 3:
#            uw_image.convert("RGB").save(sys.argv[2])
            uw_image.save(sys.argv[2])
            print("Imagen guardada con exito.")
            return True
        else:
            print("Cantidad de parametros erronea.")

    return None

if __name__=='__main__':
    exit(main())

