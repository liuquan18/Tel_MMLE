#%%
import xarray as xr
import numpy as np
#%%
def standard_index(eof_result, standard):
    # standarize the index with the tmeporal mean and std
    eof_result = eof_result.copy()
    if standard == "first":
        print(" standardizing the index with the first 10 years ...")
        years = np.unique(eof_result["time.year"])
        years = sorted(years)[:10]
        ref = eof_result["pc"].sel(time = eof_result["time.year"].isin(years))
        pc_std = (eof_result["pc"] - ref.mean(dim=("time", "ens"))) / ref.std(
            dim=("time", "ens")
        )
        eof_result["pc"] = pc_std

    # standarize the index with the temporal and ensemble mean and std
    elif standard == "temporal_ens":
        print(" standardizing the index with temporal and ensemble mean and std ...")
        eof_result = eof_result.copy()
        ref = eof_result["pc"]
        pc_std = (eof_result["pc"] - ref.mean(dim=("time", "ens"))) / ref.std(
            dim=("time", "ens")
        )
        eof_result["pc"] = pc_std
    return eof_result

# %%
models = ['MPI_GE','CanESM2','CESM1_CAM5','MK36','GFDL_CM3']

odir = '/work/mh0033/m300883/Tel_MMLE/data/'

for model in models:
    print(model)
    eof_result = xr.open_dataset(odir + model + '/EOF_result/'  + 'plev_50000_decade_mpi_JJA_none_eof_result.nc')
    eof_result = standard_index(eof_result, standard = 'first')
    eof_result.to_netcdf(odir + model + '/EOF_result/'  + 'plev_50000_decade_mpi_first_JJA_eof_result.nc')
# %%
