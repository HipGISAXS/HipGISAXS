
import os, sys, argparse
import numpy as np
import matplotlib.pyplot as plt

#from xrayutilities import io as xio


def abort(msg):
  print msg
  sys.exit(1)

def load_matrix(fname):
  ext = fname.split('.')[-1]
  if ext == 'out': return np.loadtxt(fname)
  elif ext == 'edf' or ext == 'EDF':
    e = xrayutilities.io.edf.EDFFile(fname)
    data = e.data
    return np.array(data)
  else: abort('unknown file format '+ext)

def norm(data):
  sum = 0.
  for i in range(0, len(data)):
    for j in range(0, len(data[0])):
      sum += data[i][j] * data[i][j]
  return np.sqrt(sum)

def normalize(data):
  normval = norm(data)
  for i in range(0, len(data)):
    for j in range(0, len(data[0])):
      data[i][j] = data[i][j] / normval
  return data


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



data1 = load_matrix(infile1)
data2 = load_matrix(infile2)
if maskfile is not None:
  mask = load_matrix(maskfile)
  data1 = data1 * mask
  data2 = data2 * mask

data1 = normalize(data1)
data2 = normalize(data2)

print len(data1[0]), len(data1)
print len(data2[0]), len(data2)
if len(data1[0]) != len(data2[0]) or len(data1) != len(data2):
  abort('error: the input patterns are not of the same size')

diff = np.abs(data1 - data2)
np.savetxt(outfile, diff)
