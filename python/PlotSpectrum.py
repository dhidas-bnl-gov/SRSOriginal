#!/usr/bin/env python

import sys
import matplotlib.pyplot as plt
import numpy as np


def PlotSpectrum (ifiles):

  for ifile in ifiles:
    with open(ifile, 'r') as f:

      x = []
      y = []

      for line in f:
        v = line.split()
        if (len(v) > 0 and v[0].startswith('#')):
          continue
  
        if (len(v) < 2) :
          continue
  
        x.append(float(v[0]))
        y.append(float(v[1]))

      plt.plot(x, y)

  plt.xlabel('Photon Energy [$eV$]')
  plt.ylabel('$\gamma / mm^2 / 0.1\%bw / s$')
  plt.title('Photon Spectrum')
  plt.yscale('log')
#  plt.xscale('log')

  plt.show()




if __name__ == "__main__":
  if len(sys.argv) < 2:
    print 'Usage: ', sys.argv[1], ' [ifile]s'
    exit(0)

  PlotSpectrum(sys.argv[1:])
