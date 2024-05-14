from mpi4py import MPI 
import sys
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
name = MPI.Get_processor_name()


if rank == 0:
    data = [(x+1)**x for x in range(size -1)]
    # add None at the end of the list
    data.append(None)
    print('will be scattered', data)
else:
    data = []

data = comm.scatter(data, root = 0)
# print (f'rank {rank} has data {data}' )

listdata = [data]*5

print(f'rank {rank} has list data {listdata}')

newData = comm.gather(listdata, root=0)
# flatten the list newData
# newData = [item for sublist in newData for item in sublist]

if rank == 0:
    print('master:', newData)

    # exclude any None in the elelment of newData
    flattened = [item for sublist in newData if sublist is not None for item in sublist if item is not None]
    print(f"flattened{flattened}")