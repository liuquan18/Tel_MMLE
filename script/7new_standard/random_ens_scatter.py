# %%
import src.MMLE_TEL.index_stats as index_stats
import numpy as np
import xarray as xr
import importlib
import matplotlib.pyplot as plt
import src.plots.extrc_tsurf_scatter as extrc_tsurf

# %%
# %%
def read_data(ens_size):
    eof_dir = f"/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_random/EOF_result/plev_50000_decade_temporal_ens{str(ens_size)}_eof_result.nc"
    eof_result = xr.open_dataset(eof_dir)

    tsurf_dir = "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_random/ts_processed/tsurf_mean.nc"
    tsurf = xr.open_dataset(tsurf_dir).tsurf

    return eof_result, tsurf
# %%
# extreme event count vs. tsurf
def random_scatter(ci = 'AR1',ens_size = 20):
    eof_result,tsurf = read_data(ens_size)
    print("ploting the extreme event count vs. tsurf")
    try:
        tsurf_mean = tsurf.mean(dim="ens").squeeze()
    except ValueError:
        tsurf_mean = tsurf
    tsurf_increase = tsurf_mean - tsurf_mean[0]

    ext_counts, t_surf_mean = extrc_tsurf.decadal_extrc_tsurf(
        eof_result.pc, tsurf_increase,ci = ci
    )
    extrc_tsurf_scatter = extrc_tsurf.extCount_tsurf_scatter(
        ext_counts, t_surf_mean, 
    )
    plt.savefig(f"/work/mh0033/m300883/Tel_MMLE/docs/source/plots/random_ens/ens_{str(ens_size)}" + f"_{ci}_extreme_count_tsurf.png", dpi=300)

# %%
random_scatter(ens_size=20)
# %%
random_scatter(ens_size=50)

# %%
# %%
random_scatter(ens_size=30)
# %%
random_scatter(ens_size=40)
# %%
