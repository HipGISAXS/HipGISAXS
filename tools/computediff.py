
import os, sys, argparse
import numpy as np
import matplotlib.pyplot as plt


def abort(msg):
    print msg
    sys.exit(1)

ap = argparse.ArgumentParser()
ap.add_argument('-i')   ## input matrix 1
ap.add_argument('-j')   ## input matrix 2
ap.add_argument('-o')   ## output matrix file
ap.add_argument('-m')   ## mask file
args = ap.parse_args()

infile1, infile2, outfile, maskfile = None, None, None, None
if args.i: infile1 = args.i
else: abort('error: an input data file 1 is necessary')
if args.i: infile2 = args.j
else: abort('error: an input data file 2 is necessary')
if args.o: outfile = args.o
else: abort('error: an output file name is necessary')
if args.m: maskfile = args.m


data1 = np.loadtxt(infile1)
data2 = np.loadtxt(infile2)
if maskfile is not None:
  mask = np.loadtxt(maskfile)
  data1 = data1 * mask
  data2 = data2 * mask

print len(data1[0]), len(data1)
print len(data2[0]), len(data2)
if len(data1[0]) != len(data2[0]) or len(data1) != len(data2):
  abort('error: the input patterns are not of the same size')


diff = np.abs(data1 - data2)
np.savetxt(outfile, diff)