from par_utils import *
from mpi4py import MPI

comm = MPI.COMM_WORLD

rootprint('MPI Parallel Python program. There are {:d} MPI processes.'.format(comm.size))

comm.Barrier()

allprint('Hello world from rank {:d}'.format(comm.rank))
