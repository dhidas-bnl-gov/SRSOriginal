#!/usr/bin/env python

# Import sys for command line arguments
import sys

# Correct paths to point to SRS.so library
sys.path.append('/Users/dhidas/SRS/lib')
sys.path.append('/Users/dhidas/SRS/python')
sys.path.append('/Users/dhidas/SRS/build/lib.macosx-10.11-intel-2.7')

# Import the SRS module
import SRS

# Check number of arguments is correct
if len(sys.argv) != 3:
    raise

# This section number and what to prefix output files with
SECTION_NUMBER = int(sys.argv[1])
OUTPUT_PREFIX  = sys.argv[2]

# Flux output file name
outfile = OUTPUT_PREFIX + '.' + str(SECTION_NUMBER) + '.dat'

# Create an SRS object
srs = SRS.SRS()

# Print for record
print 'Section number:', sys.argv[1]
print 'Output file:', sys.argv[2]

# Set the seed based on section number (or leave off for random device seeding)
srs.set_seed(12345 + SECTION_NUMBER)

# Add field of simple undulator or read in a file, etc.
srs.add_undulator(bfield=[0, -1, 0], period=[0, 0, 0.049], nperiods=41)

# Set the particle beam
srs.set_particle_beam(type='electron', name='beam_1', x0=[0, 0, -2], direction=[0, 0, 1], energy_GeV=3, current=0.500)

# Set the start and stop times for the calculation
srs.set_ctstartstop(0, 4)

# If you want increased precision change npoints_trajectory
#srs.set_npoints_trajectory(20001)

# Do the calculation of interest
srs.calculate_flux_rectangle(plane='XY', width=[0.01, 0.01], npoints=[11, 11], energy_eV=152., translation=[0, 0, 30], ofile=outfile, nparticles=1)

print 'done.'
