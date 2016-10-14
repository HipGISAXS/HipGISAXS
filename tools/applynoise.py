##
 #  Project:
 #
 #  File: applynoise.py
 #  Created: Sep 09, 2016
 #
 #  Author: Abhinav Sarje <asarje@lbl.gov>
 ##

import os, sys, argparse
import numpy as np


SCALE_MIN = 0
SCALE_MAX = 2**20

def abort(msg):
  print msg
  sys.exit(1)

ap = argparse.ArgumentParser()
ap.add_argument('-i')
ap.add_argument('-o')
ap.add_argument('-n')
args = ap.parse_args()

infile, outfile = None, None
if args.i: infile = args.i
else: abort('error: an input data file is necessary')
if args.o: outfile = args.o
else: abort('error: an output file name is necessary')
if args.n: N = int(args.n)

## load the data
data = np.loadtxt(infile)

## find the min and max values (nonzero)
#xmin = np.min(data[np.nonzero(data)])
xmin = 0.
xmax = np.max(data[np.nonzero(data)])

## scale to [0, 2^20] unsigned integers (20-bit)
ymin, ymax = SCALE_MIN, SCALE_MAX
print '++ scaling from [', xmin, xmax, '] to [', ymin, ymax, ']'
scaled_data = np.array([[((ymax - ymin) * (x - xmin) / (xmax - xmin) + ymin) for x in xrow] for xrow in data])
np.savetxt('scaled_'+infile, scaled_data, fmt='%d')

print '++ applying poisson noise'
#noisy_data = np.random.poisson(scaled_data)
noisy_data = np.array([ [ np.random.poisson(x) for x in xrow ] for xrow in scaled_data ])
np.savetxt(outfile, noisy_data, fmt='%d')

