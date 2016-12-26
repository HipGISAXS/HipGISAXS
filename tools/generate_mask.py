##
 #  Project:
 #
 #  File: generate_mask.py
 #  Created: Oct 31, 2016
 #
 #  Author: Abhinav Sarje <asarje@lbl.gov>
 ##


import sys, argparse
import numpy as np


XSIZE = 327
YSIZE = 332

MASK_X_BEGIN = 145
MASK_X_END = 182

MASK_Y_BEGIN = 0
MASK_Y_END = 332

MASK_BEGINS = [ (145, 0), (0, 270) ]
MASK_ENDS = [ (182, 332), (327, 332) ]


def abort(msg):
  print msg
  sys.exit(1)


inmask = 'mask.out'
outmask = 'new_mask.out'

mask_data = np.loadtxt(inmask)
if XSIZE != len(mask_data[0]) or YSIZE != len(mask_data): abort('mask size not correctly set')

for j in range(0,YSIZE):
  for i in range(0,XSIZE):
    for ( x_beg, y_beg ), ( x_end, y_end ) in zip(MASK_BEGINS, MASK_ENDS):
      if i >= x_beg and i <= x_end and j >= y_beg and j <= y_end: mask_data[j][i] = 0

np.savetxt(outmask, mask_data, fmt='%d')
