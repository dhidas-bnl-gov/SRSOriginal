#!/usr/bin/env python

import sys
import numpy
import matplotlib.pyplot as plt

# The paths below are just convenient for me on my mac for testing...
sys.path.append('./')
sys.path.append('./lib')
sys.path.append('./build/lib.macosx-10.11-intel-2.7')

# Import the SRS module
import SRSP

# Create a new SRS object
srsp0 = SRSP.SRSP()

# Fille some variables
for i in range(10):
  srsp0.set_ctstart(i)
  #srsp0.ctstart = i
  print srsp0.ctstart

  srsp0.set_ctstop(i+10)
  print srsp0.ctstop
  #print srsp0.get_ctstop()


print ''
print ''
print ''


for i in range(10):
  srsp0.set_ctstartstop(i, i + 10)
  print srsp0.get_ctstart(), srsp0.get_ctstop()

srsp0.set_npoints_trajectory(123)
print srsp0.get_npoints_trajectory()

srsp0.add_magnetic_field('epu_linear.dat', 'ZBxByBz')

Z = []
By = []
for z in numpy.linspace(-1, 1, 10000):
  Z.append(z)
  By.append(srsp0.get_bfield([0, 0, z])[1])


plt.plot(Z, By)
plt.xlabel('Z [m]')
plt.ylabel('$B_Y$ [T]')
plt.title('Magnetic Field')
plt.show()
