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
def read_data( ens_size,model = 'MPI_GE_onepct',plev = 50000,standard = 'first',tsurf = 'ens_fld_year_mean'):
    eof_dir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}_random/EOF_result/plev_{plev}_decade{str(ens_size)}_{standard}_eof_result.nc"
    eof_result = xr.open_dataset(eof_dir)

    tsurf_dir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}_random/ts_processed/{tsurf}.nc"
    tsurf = xr.open_dataset(tsurf_dir).tsurf

    extrc_dir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}_random/extreme_count/plev_{plev}_decade{str(ens_size)}_{standard}_extre_counts.nc"
    try:
        extrc = xr.open_dataset(extrc_dir).pc
    except FileNotFoundError:
        extrc = None
    return eof_result, tsurf, extrc

# %%
# count the extreme event and tsurf
def extreme_counts (eof_result, tsurf, extrc = None):

    tsurf_increase = tsurf - tsurf[0]
    if extrc is None:
        ext_counts, t_surf_mean = extrc_tsurf.decadal_extrc_tsurf(
            eof_result.pc, temp = tsurf_increase,ci = 'bootstrap'
        )
    else:
        ext_counts, t_surf_mean = extrc_tsurf.decadal_extrc_tsurf(
            eof_result.pc, ext_counts_xr = extrc, temp=tsurf_increase, ci = 'bootstrap'
        )
    return ext_counts, t_surf_mean

# %%
from multiprocessing import Pool

ens_sizes = np.arange(20,101,10)
tsurfs = ["ens_fld_year_mean","NA_tsurf","tropical_arctic_gradient"]

def process_data(ens_size, tsurf):
    eof_result, temperature, extrc= read_data(ens_size,tsurf = tsurf)
    # close the files red
    temperature.close()
    
    extr_counts, tsurf_mean = extreme_counts(eof_result, temperature, extrc=extrc)
    try:
        extr_counts.to_netcdf(f"/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_random/extreme_count/plev_50000_first_extreme_counts_{str(ens_size)}_{tsurf}.nc")
    except PermissionError:
        pass
    tsurf_mean.to_netcdf(f"/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_random/extreme_count/{tsurf}.nc")

if __name__ == '__main__':
    with Pool() as p:
        p.starmap(process_data, [(ens_size, tsurf) for ens_size in ens_sizes for tsurf in tsurfs])
    


# %%
