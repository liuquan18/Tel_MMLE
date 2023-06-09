#%%
import xarray as xr
import numpy as np
import src.extreme.extreme_ci as extreme
import matplotlib.pyplot as plt
from pathlib import Path
# %%
import importlib
importlib.reload(extreme)

# %%
def standard_by_first(eof):
    pcs = eof.pc
    pcs_std = pcs.isel(time = slice(0,10)).std(dim = ("time","ens"))
    pcs_mean = pcs.isel(time = slice(0,10)).mean(dim = ("time","ens"))
    std_pcs = (pcs - pcs_mean)/pcs_std
    eof['pc'] = std_pcs
    return eof

#%%
def standard_models(model):
    odir = "/work/mh0033/m300883/Tel_MMLE/data/"
    EOF_dir = Path(odir) / model / "EOF_result/"
    for eof_result_name in EOF_dir.glob("*_none_eof_result.nc"):
        new_name = eof_result_name.name.replace("none","first")
        # check if the new file exists
        if (EOF_dir / new_name).exists():
            continue
        eof = xr.open_dataset(eof_result_name)
        std_eof = standard_by_first(eof)
        std_eof.to_netcdf(EOF_dir / new_name)
        print(eof_result_name.name + " is done")

#%%

import multiprocessing

def process_model(model):
    standard_models(model)

models = ["CanESM2","CESM1_CAM5","GFDL_CM3","MK36","MPI_GE","MPI_GE_onepct","MPI_GE_onepct_random"]

# Create a pool of worker processes
pool = multiprocessing.Pool()

# Process each model in parallel
pool.map(process_model, models)

# Close the pool of worker processes
pool.close()
pool.join()



#%%
if __name__ == "__main__":
    # %%
    trop_eof = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct/EOF_result/troposphere_ind_decade_none_eof_result.nc")
    std_pcs = standard_by_first(trop_eof).pc
    # %%
    tsurf = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct/ts_processed/ens_fld_year_mean.nc").tsurf
    # %%
    tsurf = tsurf.mean(dim = "ens").squeeze()
    tsurf_increase = tsurf - tsurf.isel(time = 0)

    # %%
    import src.plots.extrc_tsurf_scatter as extrc_tsurf

    # %%
    ext_counts, t_surf_mean = extrc_tsurf.decadal_extrc_tsurf(
        std_pcs, tsurf_increase, ci='bootstrap'
    )
    # %%
    extrc_tsurf_scatter = extrc_tsurf.extCount_tsurf_scatter(
        ext_counts.sel(plev = 50000).squeeze(), t_surf_mean
    )
    # %%
    # profile
    first_pc = std_pcs.isel(time=slice(0, 10))
    last_pc = std_pcs.isel(time=slice(-10, None))
    # %%
    print("calculating the extreme event count")
    # extreme counts of first and last 10 decades
    first_count = extreme.extreme_count_xr(first_pc, ci='bootstrap')
    last_count = extreme.extreme_count_xr(last_pc, ci='bootstrap')
    # %%
    extreme_profile = extreme.extreme_count_profile(
        first_count, last_count, colored=False,xlim = (45,130)
    )
    # %%
