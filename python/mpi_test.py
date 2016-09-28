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

srs.add_particle_beam(type='electron', name='beam_0', energy_GeV=3, direction=[0, 0, 1], x0=[0, 0, -1], current=0.5)
srs.set_ctstartstop(0, 2)
srs.add_undulator(bfield=[0, 1, 0], period=[0, 0, 0.050], nperiods=33)


npoints = 500


if rank == 0:
    # Just sit here and collect data

    # Weight for each spectrum in summation
    weight = 1. / size

    # Collect data from other processes
    for i in range(1, size):
        data = comm.recv(source=ANY_SOURCE, tag=11)

        # Sum the spectra (internally this is a compensated sum)
        srs.add_to_spectrum(spectrum=data, weight=weight)

    plot_spectrum(srs.get_spectrum())

else:
    # all other process send their result

    data = srs.calculate_spectrum(obs=[0, 0, 30], energy_range_eV=[100, 1000], npoints=npoints, nthreads=1)
    comm.send(data, dest=0, tag=11)

