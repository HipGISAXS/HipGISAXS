
import os, sys, argparse
import numpy as np
import matplotlib.pyplot as plt


def abort(msg):
	print msg
	sys.exit(1)

def norm_l2(data, mask):
	datanorm = 0.
	d2 = data * data
	if mask is not None: d2 = d2 * mask
	datanorm = d2.sum()
	return np.sqrt(datanorm)

def dot_prod(data1, data2, mask):
	d2 = data1 * data2
	if mask is not None: d2 = d2 * mask
	dotprod = d2.sum()
	return dotprod

ap = argparse.ArgumentParser()
ap.add_argument('-i', help='input file')
ap.add_argument('-o', help='output file')
ap.add_argument('-m', help='mask file')
ap.add_argument('-r', help='reference data file')
ap.add_argument('-axes', help='plot axes labels?')
ap.add_argument('-cbar', help='plot color bar?')
ap.add_argument('-log', help='plot on log scale?')
args = ap.parse_args()

infile, outfile, maskfile, reffile = None, None, None, None
axes, log, cbar = True, True, True
if args.i: infile = args.i
else: abort('error: an input data file is necessary')
if args.o: outfile = args.o
else: abort('error: an output file name is necessary')
if args.m: maskfile = args.m
if args.r: reffile = args.r
if args.axes: axes = (int(args.axes) != 0)
if args.cbar: cbar = (int(args.cbar) != 0)
if args.log: log = (int(args.log) != 0)

plotfile = outfile

data = np.loadtxt(infile)
if maskfile is not None:
  mask = np.loadtxt(maskfile)
  data = data * mask
else: mask = None
	
print len(data[0]), len(data)


## normalize with the reference data (if given)
if reffile is not None:
	rdata = np.loadtxt(reffile)
	if mask is not None: rdata = rdata * mask
	dotprod = dot_prod(data, rdata, mask)
	datnorm = norm_l2(data, mask)
	data = data * dotprod / (datnorm * datnorm)
	# MAXVAL = 2 ** 20
	# dmin, dmax = data.min(), data.max()
	# data = MAXVAL * (data - dmin) / (dmax - dmin)
	# data = np.array([ map(int, d) for d in data ])

# none nearest bilinear bicubic spline16 spline36 hanning hamming hermite kaiser quadric catrom gaussian bessel mitchell sinc lanczos

if log: im = plt.imshow(np.log(data[:,:]), interpolation='gaussian', cmap='RdYlBu_r', origin='upper')
else: im = plt.imshow(data[:,:], interpolation='gaussian', cmap='RdYlBu_r', origin='upper')

if axes:
	plt.xlabel('qy (nm$^{-1}$)')
	plt.ylabel('qz (nm$^{-1}$)')
if cbar: plt.colorbar()

plt.gcf().set_size_inches(10, 5)
plt.savefig(plotfile, bbox_inches='tight')