#!/usr/bin/env python

import sys
import matplotlib.pyplot as plt
import numpy as np


def PlotFlux (fileName) :
  f = open(fileName, 'r')

  x = []
  y = []
  z = []

  with open(fileName, 'r') as f:
    for line in f:
      v = line.split()
      if (len(v) > 0 and v[0].startswith('#')):
        continue
  
      if (len(v) < 3) :
        continue
  
      x.append(float(v[0]))
      y.append(float(v[1]))
      z.append(float(v[2]))


  nx = len(np.unique(x))
  ny = len(np.unique(y))

  plt.hist2d(x, y, bins=[nx, ny], weights=z)
  plt.colorbar()

  plt.xlabel('Axis 1 [$m$]')
  plt.ylabel('Axis 2 [$m$]')
  plt.title('Flux ' +'[$\gamma / mm^2 / 0.1bw / s$]')

  #plt.axes().set_aspect('equal', 'datalim')
  #plt.set_cmap('gray') # 'hot' is also nice

  plt.show()




if __name__ == "__main__":
  if len(sys.argv) is not 2:
    print 'Usage: ', sys.argv[1], ' [InFile]'
    exit(0)

  PlotFlux(sys.argv[1])
