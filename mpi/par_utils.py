from mpi4py import MPI

def rootprint(str="", end="\n", comm=MPI.COMM_WORLD):
    if comm.rank == 0:
        print(str+end)


def allprint(str="", end="\n", comm=MPI.COMM_WORLD):
    print('Rank ',comm.rank, str+end)
