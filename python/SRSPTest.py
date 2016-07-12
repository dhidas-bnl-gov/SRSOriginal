#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
from math import sin, cos

# The paths below are just convenient for me on my mac for testing...
sys.path.append('./')
sys.path.append('./lib')
sys.path.append('./build/lib.macosx-10.11-intel-2.7')

# Import the SRS module
import SRS

# Create a new SRS object
srs0 = SRS.SRS()

print srs0.pi()
exit(0)


def myfunc (X):
  z = X[2]
  return [sin(5*z), cos(5*z), 0]



srs0.add_magnetic_field_function(myfunc)



if False:
  Z = np.linspace(-4, 4, 2000)
  B = [srs0.get_bfield([0, 0, z])[0] for z in Z]
  # Bx = [item[0] for item in B]
  plt.plot(Z, B)
  plt.xlabel('Z [m]')
  plt.ylabel('Bx [T]')
  plt.title('Magnetic Field')
  plt.show()







srs0.set_ctstartstop(-2, 2)
srs0.set_npoints_trajectory(20001)

#srs0.add_magnetic_field('epu_linear.dat', 'ZBxByBz')
srs0.add_particle_beam('electron', 'beam_0', [0, 0, 0], [0, 0, 1], 3, -2, 0.500, 0.1)
#srs0.add_particle_beam('electron', 'beam_1', [0, 0, -2], [0, 0, 1], 3, -2, 0.500, 0.9)

srs0.set_new_particle()


if False:
  for i in xrange(10000):
    srs0.set_new_particle()



if False:
  srs0.calculate_trajectory()
  trajectory = srs0.get_trajectory()

  Z = [item[0][2] for item in trajectory]
  X = [item[0][0] for item in trajectory]
  plt.plot(Z, X)
  plt.xlabel('Z [m]')
  plt.ylabel('X [m]')
  plt.title('Particle Trajectory')
  plt.show()

power_density = srs0.calculate_power_density_rectangle2(width_x1=0.008, nx1=51, width_x2=0.004, nx2=51, rotations=[0, 0, 30], translation=[0, 0, 0], x0x1x2=[[1,1,1],[2,2,2],[3,3,3]], plane="XY")
#power_density = srs0.calculate_power_density_rectangle2(width_x1=0.008, nx1=51, width_x2=0.004, nx2=51, rotations=[0, 0, 30], translation=[0, 0, 0], points=[[1,1,1],[2,2,2],[3,3,3]], plane="XY")

exit(0)




if False:
  #srs0.calculate_spectrum([0, 0, 30], 100, 500, 1000)
  srs0.calculate_spectrum_from_list([0, 0, 30], range(100, 500))
  Spectrum = srs0.get_spectrum()

  X = [item[0] for item in Spectrum]
  Y = [item[1] for item in Spectrum]
  plt.plot(X, Y)
  plt.xlabel('Energy [eV]')
  plt.ylabel('Photons [$\gamma / mm^2 / 0.1bw / s$]')
  plt.title('Spectrum')
  plt.show()


if False:
  power_density = srs0.calculate_power_density_rectangle('XY', 0.008, 51, 0.004, 51, [0, 0, 30], 1)
  X = [item[0][0] for item in power_density]
  Y = [item[0][1] for item in power_density]
  P = [item[1]    for item in power_density]

  NX = len(np.unique(X))
  NY = len(np.unique(Y))

  plt.hist2d(X, Y, bins=[NX, NY], weights=P)
  plt.colorbar()

  plt.xlabel('Axis 1 [$m$]')
  plt.ylabel('Axis 2 [$m$]')
  plt.title('Power Density ' +'[$W / mm^2$]')

  plt.show()


if True:
  flux = srs0.calculate_flux_rectangle(135, 'XY', 0.008, 51, 0.004, 51, [0, 0, 30], 1)
  X = [item[0][0] for item in flux]
  Y = [item[0][1] for item in flux]
  P = [item[1]    for item in flux]

  NX = len(np.unique(X))
  NY = len(np.unique(Y))

  plt.hist2d(X, Y, bins=[NX, NY], weights=P)
  plt.colorbar()

  plt.xlabel('Axis 1 [$m$]')
  plt.ylabel('Axis 2 [$m$]')
  plt.title('Spectral Flux' +'[$W / mm^2$]')

  plt.show()





