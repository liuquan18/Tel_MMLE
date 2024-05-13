# %% === Import ===
import numpy as np
import xarray as xr
import mpi4py.MPI as MPI
import glob
from src.Teleconnection.spatial_pattern import project_field
from src.MMLE_TEL.index_generator import read_data
from src.MMLE_TEL.index_generator import standard_index
import pandas as pd
import sys
import logging
#%%
# logging configure
logging.basicConfig(level=logging.INFO)
# %%
num = int(sys.argv[1]) # get the index for the node
model = sys.argv[2] # get the model name
plev = 50000
logging.info(f"Node {num} is working on {model}")

# === mpi4py ===
try:
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()  # 0,1,2,3,4,5,6,7,8,9
    npro = comm.Get_size()  # 10
except:
    print("::: Warning: Proceeding without mpi4py! :::")
    rank = 0
    npro = 1


#%%
# define the path
data_dir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/zg_JJA/"
index_dir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/EOF_result/"

#%%
# list all the files under data_dir
files = glob.glob(data_dir + "*.nc")
files.sort()

# read eof data from the whole period
eof_path = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/EOF_result/plev_50000_all_first_JJA_eof_result.nc"
eof_result = xr.open_dataset(eof_path)
eof = eof_result.eof.sel(mode = 'NAO').squeeze()

files_single_core = np.array_split(files, npro)[rank]

index_allens = []

for step, file in enumerate(files_single_core):
    filename = file.split("/")[-1]
    logging.info(f"Processing {filename} on node {num} step {step}/{len(files_single_core)}")
    # read zg_data
    data = xr.open_dataset(file)
    # select plev
    zg_data = data.sel(plev = plev)
    try:
        zg_data = zg_data.var156
    except AttributeError:
        zg_data = zg_data.zg
    
    # time to datetime
    try:
        zg_data["time"] = zg_data.indexes["time"].to_datetimeindex()
    except AttributeError:
        zg_data["time"] = pd.to_datetime(zg_data.time)
    zg_data = zg_data.squeeze()

    # project the field

    index = project_field(zg_data, eof, dim = 'time', standard=False)

    index_allens.append(index)

logging.info ("collecting the index from all processes")
# gather the index from all processes
index_allens = comm.gather(index_allens, root=0)

logging.info ("index start to be standardized")
# standardize the index
if rank == 0:

    def standard_index(index_allens):
        # standardize the index with the mean and std of the first 10 years, all ens
        ref_mean = index_allens.isel(time = slice(0,10)).mean(dim = ('time','ens'))
        ref_std = index_allens.isel(time = slice(0,10)).std(dim = ('time','ens'))

        standarded = (index_allens - ref_mean) / ref_std
        return standarded
    

    index_allens = [item for sublist in index_allens for item in sublist]
    index_allens = xr.concat(index_allens, dim = 'ens')
    index_allens = standard_index(index_allens)

    # save the index
    index_allens.to_netcdf(index_dir + "plev_50000_all_first_JJA_projected_index.nc")
    logging.info ("index saved")

else:
    pass