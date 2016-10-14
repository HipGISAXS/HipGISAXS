
import os, sys, argparse
import numpy as np
import matplotlib.pyplot as plt


def abort(msg):
    print msg
    sys.exit(1)

ap = argparse.ArgumentParser()
ap.add_argument('-i')
ap.add_argument('-o')
ap.add_argument('-m')
args = ap.parse_args()

infile, outfile, maskfile = None, None, None
if args.i: infile = args.i
else: abort('error: an input data file is necessary')
if args.o: outfile = args.o
else: abort('error: an output file name is necessary')
if args.m: maskfile = args.m

plotfile = outfile

data = np.loadtxt(infile)
if maskfile is not None:
    mask = np.loadtxt(maskfile)
    data = data * mask
print len(data[0]), len(data)

# none nearest bilinear bicubic spline16 spline36 hanning hamming hermite kaiser quadric catrom gaussian bessel mitchell sinc lanczos

im = plt.imshow(np.log(data[:,:]), interpolation='none', cmap='RdYlBu_r', origin='upper')
#im = plt.imshow(np.log(data[:,:]), interpolation='gaussian', cmap='RdYlBu_r', origin='upper')
#im = plt.imshow(data[:,:], interpolation='gaussian', cmap='RdYlBu_r', origin='upper')

plt.xlabel('qy (nm$^{-1}$)')
plt.ylabel('qz (nm$^{-1}$)')
#plt.colorbar()
plt.gcf().set_size_inches(10, 5)
plt.savefig(plotfile, bbox_inches='tight')
