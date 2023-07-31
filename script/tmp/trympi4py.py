from mpi4py import MPI 
import sys
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
name = MPI.Get_processor_name()

shared = (rank + 1) * 5

if rank == 0:
    data = shared
    comm.send(data, dest=1, tag=11)
    print(f"rank {rank} sent {data} to rank 1")

if rank == 1:
    data = comm.recv(source=0, tag=11)
    print(f"rank {rank} received {data} from rank 0")