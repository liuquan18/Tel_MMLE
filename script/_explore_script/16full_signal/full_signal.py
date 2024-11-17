#%%
import xarray as xr
import numpy as np
import glob

from src.Teleconnection.spatial_pattern import project_field
from src.MMLE_TEL.index_generator import read_data
import matplotlib.pyplot as plt
#%%

from src.compute.slurm_cluster import init_dask_slurm_cluster

#%%

client, cluster = init_dask_slurm_cluster(scale = 10, processes = 2, walltime="05:30:00", memory="200GB")


# %%
model = 'MPI_GE'
#%%
def read_zg_data(model, plev = 50000):
        # read gph data
        odir = odir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/"
        data_JJA = []
        for month in ["Jun", "Jul", "Aug"]:
            print(f"reading the gph data of {month} ...")
            zg_path = odir + "zg_" + month + "/"
            data_JJA.append(read_data(zg_path, plev=plev, remove_ens_mean=False))
        data = xr.concat(data_JJA, dim="time").sortby("time")

        return data

# %%
data = read_zg_data(model)
# %%

eof_first = xr.open_dataset(f"/work/mh0033/m300883/Tel_MMLE/data/{model}/EOF_result/first_pattern_projected.nc").__xarray_dataarray_variable__
eof_last = xr.open_dataset(f"/work/mh0033/m300883/Tel_MMLE/data/{model}/EOF_result/last_pattern_projected.nc").__xarray_dataarray_variable__

#%%
eof_first.load()
eof_last.load()

# %%
data_first = data.sel(time = slice('1850-06-01', '1859-08-31')) # time may be different
data_last = data.sel(time = slice('2090-06-01', '2099-08-31'))
#%%
data_first = data_first.stack(com = ('time','ens'))
data_last = data_last.stack(com = ('time','ens'))
# %%
ppc_first = project_field(data_first, eof_first)
# %%
ppc_last = project_field(data_last, eof_first)
# %%
fig, ax = plt.subplots()
ppc_first.plot.hist(ax = ax)
ppc_last.plot.hist(ax = ax)
plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/mechism/NAO_full_signal_{model}.png")

# %%
ppc_first.name = 'pc'
ppc_last.name = 'pc'

ppc_first.to_netcdf(f"/work/mh0033/m300883/Tel_MMLE/data/{model}/full_signal/first_pattern_projected.nc")
ppc_last.to_netcdf(f"/work/mh0033/m300883/Tel_MMLE/data/{model}/full_signal/last_pattern_projected.nc")
# %%
first_ano = (ppc_first - ppc_first.mean())/ppc_first.std()
last_ano = (ppc_last - ppc_first.mean())/ppc_first.std()

#%%
first_ano.to_netcdf(f"/work/mh0033/m300883/Tel_MMLE/data/{model}/full_signal/first_pattern_projected_ano.nc")
last_ano.to_netcdf(f"/work/mh0033/m300883/Tel_MMLE/data/{model}/full_signal/last_pattern_projected_ano.nc")

# %%
fig, ax = plt.subplots()
first_ano.plot.hist(ax = ax, alpha = 0.7,bins=np.arange(-4, 4.1, 0.5),color = 'k', label = 'first10')
last_ano.plot.hist(ax = ax, alpha = 0.7, bins=np.arange(-4, 4.1, 0.5),color = 'r', label = 'last10')

# vline at mean
plt.axvline(first_ano.mean(), color='k', linestyle='dashed', linewidth=1.5)
plt.axvline(last_ano.mean(), color='r', linestyle='dashed', linewidth=1.5)

# title
plt.title("full change in NAO")

plt.legend()
# plt.savefig(f"/work/mh0033/m300883/Tel_MMLE/docs/source/plots/mechism/NAO_full_signal_ano_{model}.png")

# %%
