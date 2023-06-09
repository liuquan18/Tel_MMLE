# %%
import src.MMLE_TEL.index_stats as index_stats
import numpy as np
import importlib
import matplotlib.pyplot as plt
import src.plots.extrc_tsurf_scatter as extrc_tsurf
import xarray as xr

#%%
import importlib

importlib.reload(index_stats)


# %%
def read_data( ens_size,model = 'MPI_GE_onepct',plev = 50000,standard = 'first'):
    eof_dir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}_random/EOF_result/plev_{plev}_decade{str(ens_size)}_{standard}_eof_result.nc"
    eof_result = xr.open_dataset(eof_dir)

    tsurf_dir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}_random/ts_processed/ens_fld_year_mean.nc"
    tsurf = xr.open_dataset(tsurf_dir).tsurf
    tsurf = tsurf.mean(dim = 'ens').squeeze()

    return eof_result, tsurf

# %%
# count the extreme event and tsurf
def extreme_counts (eof_result, tsurf):
    tsurf_increase = tsurf - tsurf[0]
    ext_counts, t_surf_mean = extrc_tsurf.decadal_extrc_tsurf(
        eof_result.pc, tsurf_increase,ci = 'bootstrap'
    )

    ds = xr.Dataset(
        {
            "extreme_counts": ext_counts,
            "tsurf": t_surf_mean,
        }
    )
    return ds

# %%
from multiprocessing import Pool

ens_sizes = np.arange(20,101,10)

def process_ens_size(ens_size,standard = 'first'):
    plev = 50000
    eof_result, tsurf = read_data(ens_size,plev = plev)
    ds = extreme_counts(eof_result, tsurf)
    ds.to_netcdf(f"/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_random/extreme_count/plev_{plev}_{standard}_extreme_counts_tsurf_{str(ens_size)}.nc")

with Pool() as p:
    p.map(process_ens_size, ens_sizes)

# %%
