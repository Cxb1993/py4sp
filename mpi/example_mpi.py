from par_utils import *
from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD

rootprint('MPI Parallel Python program. There are {:d} MPI processes.'.format(comm.size))
comm.Barrier()

allprint('Hello world from rank {:d}'.format(comm.rank))
comm.Barrier()

# Test some collective comms
random_number = np.random.rand(1)
allprint('My random number is {:4.4f}'.format(random_number[0]))
comm.Barrier()
# broadcast
comm.Bcast([random_number, MPI.DOUBLE])

comm.Barrier()
rootprint('--------------------------------------------------')
allprint('My random number is {:4.4f}'.format(random_number[0]))
comm.Barrier()

# reduce
#
## Test some point-to-point comms
