# How to run this code:
#   mpiexec -n 5 python python/mpi_test.py

import sys
sys.path.append('/Users/dhidas/SRS/lib')
sys.path.append('/Users/dhidas/SRS/python')

# Import the SRS and helper modules
import SRS
from SRS_BasicPlots import *

# MPI imports
from mpi4py import MPI
from mpi4py.MPI import ANY_SOURCE
import numpy

# Common communication, rank, size
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


# Get an SRS object
srs = SRS.SRS()

srs.set_particle_beam(type='electron',
                      name='beam_0',
                      x0=[0, 0, -1.5],
                      direction=[0, 0, 1],
                      energy_GeV=3,
                      current=0.500,
                      sigma_energy_GeV=0.001*3.,
                      horizontal_direction=[1, 0, 0],
                      beta=[20.1, 3.4],
                      emittance=[0.55e-9, 0.008e-9],
                      lattice_center=[0, 0, 0]
                     )


srs.set_ctstartstop(0, 2)
srs.add_undulator(bfield=[0, 1, 0], period=[0, 0, 0.050], nperiods=33)


npoints = 200


if rank == 0:
    srs.set_new_particle(particle='ideal')
    spectrum_se = srs.calculate_spectrum(obs=[0, 0, 30], energy_range_eV=[400, 460], npoints=npoints, nthreads=1)
    # Now just sit here and collect data

    # Weight for each spectrum in summation
    weight = 1. / size

    # Collect data from other processes
    srs_sum = SRS.SRS()
    for i in range(1, size):
        data = comm.recv(source=ANY_SOURCE, tag=11)

        # Sum the spectra (internally this is a compensated sum)
        srs_sum.add_to_spectrum(spectrum=data, weight=weight)

    plot_spectra([spectrum_se, srs_sum.get_spectrum()], ['single-electron', 'multi-electron'])

else:
    # all other process send their result

    data = srs.calculate_spectrum(obs=[0, 0, 30], energy_range_eV=[400, 460], npoints=npoints, nthreads=1, nparticles=20)
    comm.send(data, dest=0, tag=11)

