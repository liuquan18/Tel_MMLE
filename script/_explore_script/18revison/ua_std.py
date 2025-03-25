#%%
import xarray as xr
import src.compute.slurm_cluster as slurm_cluster

#%%
client, cluster = slurm_cluster.init_dask_slurm_cluster(scale = 2, processes = 1, memory = "200GB", walltime = "02:00:00")


# %%
odir = "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/u/"
ua = xr.open_mfdataset(odir + "*.nc", combine='nested', concat_dim='ens', parallel=True)

# change chunk size
ua = ua.chunk({'time': 30, 'lat': -1, 'lon': -1})

ua_first = ua.sel(time = slice('1850','1859'), plev = 50000).var131
ua_last = ua.sel(time = slice('2090', '2099'), plev = 50000).var131

ua_first_std = ua_first.std(dim=('time','ens'))
ua_last_std = ua_last.std(dim=('time','ens'))

ua_first_std.to_netcdf("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/zg_first_last/ua_first_std.nc")
ua_last_std.to_netcdf("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/zg_first_last/ua_last_std.nc")
