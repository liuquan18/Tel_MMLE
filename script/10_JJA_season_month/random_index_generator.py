#%%
#%%
import src.MMLE_TEL.index_generator as index_generate
import xarray as xr
import numpy as np
import sys
# %%
import importlib
importlib.reload(index_generate)


# %%
def decompose_random(ens_size,model = 'MPI_GE_onepct', fixedPattern = 'decade_mpi',plev = 50000):
    index_gen= index_generate.decompose_plev_random_ens(base_model=model,fixedPattern =fixedPattern, ens_size=ens_size,standard='first',plev=plev)
    index_gen.save_result()

#%%
ens_sizes = [20, 30, 40, 50,60,70,80,90,100]
model = 'MPI_GE_onepct'


#%%
num = int(sys.argv[1])
t1 = int(sys.argv[2])
t2 = int(sys.argv[3])

#%%
num = 1
#%%
ens_size = ens_sizes[num - 1]
print("++++++++++++++++++++++++++")
print(f"node_num {num} is running ens_size {ens_size}")
decompose_random(ens_size=ens_size,model = model)
