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
srs = SRS.SRS()



DATADIR = './Data_CSX'



srs.add_magnetic_field(DATADIR + '/mag_low_beta.lab.dat', 'ZBxByBz', [0, 0, 0], [0, 0, 0], [1, 1, -1])
#srs.add_magnetic_field(DATADIR + '/EPU1_1141a/B_G11.5_M0_Ph0_X0_Y0.txt', 'ZBxByBz', [0, 0, 0], [0, 0, -1.195], [1e-3])
#srs.add_magnetic_field(DATADIR + '/EPU2_1141b/B_G11.5_M0_Ph0_X0_Y0.txt', 'ZBxByBz', [0, 0, 0], [0, 0,  1.290], [1e-3])
srs.add_magnetic_field(DATADIR + '/EPU1_1141a/B_G11.5_M0_Ph12.3_X0_Y0.txt', 'ZBxByBz', [0, 0, 0], [0, 0, -1.195], [1e-3])
srs.add_magnetic_field(DATADIR + '/EPU2_1141b/B_G11.5_M0_Ph12.3_X0_Y0.txt', 'ZBxByBz', [0, 0, 0], [0, 0,  1.290], [1e-3])
#srs.add_magnetic_field_gaussian([0, -0.00099, 0], [0, 0, 0], [0, 0, 0.05])
#srs.add_magnetic_field_gaussian([0.00003, 0, 0], [0, 0, 0], [0, 0, 2])


# Dipole corrections for undulators
srs.add_magnetic_field_uniform([+0.00002, -0.00019, 0], [0, 0, 2.05520], [0, 0, -1.19500])
srs.add_magnetic_field_uniform([-0.00001, -0.00023, 0], [0, 0, 2.05520], [0, 0, +1.29000])

# Just to slightly correct the exit angle
srs.add_magnetic_field_gaussian([0.00035, 0, 0], [0, 0, 2.40833], [0, 0, 0.050])

if False:
    Z = np.linspace(-7.6, 3, 10000)
    B = [srs.get_bfield([0, 0, z])[1] for z in Z]
    # Bx = [item[0] for item in B]
    plt.plot(Z, B)
    plt.xlabel('Z [m]')
    plt.ylabel('Bx [T]')
    plt.title('Magnetic Field')
    plt.show()





srs.set_ctstartstop(-4, 10)
srs.set_npoints_trajectory(20001)


angle  = 0.00023
offset = 0.0005
#angle  = 0.000
#offset = 0.000

z0 = -4
epu1_center = -1.195

initial_position  = [0, angle * (z0 - epu1_center) + offset, z0]
initial_direction = [0, 1 * sin(angle), 1 * cos(angle)]

print initial_position, initial_direction

srs.add_particle_beam('electron', 'beam_0', initial_position, initial_direction, 3, 0, 0.250, 1)

srs.set_new_particle()




if False:
    srs.calculate_trajectory()
    trajectory = srs.get_trajectory()

    Z = [item[0][2] for item in trajectory]
    X = [item[0][0] for item in trajectory]
    plt.plot(Z, X)
    plt.xlabel('Z [m]')
    plt.ylabel('X [m]')
    plt.title('Particle Trajectory')
    plt.show()




if False:
    #srs.calculate_spectrum([0, 0, 30], 100, 500, 1000)
    srs.calculate_spectrum_from_list([0, 0, 30], range(100, 2000))
    Spectrum = srs.get_spectrum()

    X = [item[0] for item in Spectrum]
    Y = [item[1] for item in Spectrum]
    plt.plot(X, Y)
    plt.xlabel('Energy [eV]')
    plt.ylabel('Photons [$\gamma / mm^2 / 0.1bw / s$]')
    plt.title('Spectrum')
    plt.show()


if True:
#    power_density = srs.calculate_power_density_rectangle('XY', 0.050, 51, 0.050, 51, [0, 0, 30], 1)
    power_density = srs.calculate_power_density_rectangle('XZ', 0.080, 51, 2.40833, 51, [0, 0.004, 1.290], 1)
    X = [item[0][0] for item in power_density]
    Y = [item[0][2] for item in power_density]
    P = [item[1]    for item in power_density]

    NX = len(np.unique(X))
    NY = len(np.unique(Y))

    plt.hist2d(X, Y, bins=[NX, NY], weights=P)
    plt.colorbar()

    plt.xlabel('Axis 1 [$m$]')
    plt.ylabel('Axis 2 [$m$]')
    plt.title('Power Density ' +'[$W / mm^2$]')

#    plt.show()


if False:
    flux = srs.calculate_flux_rectangle(135, 'XY', 0.008, 51, 0.004, 51, [0, 0, 30], 1)
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





