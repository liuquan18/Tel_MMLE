# plot the extreme_count v.s tsurf for all models
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

import multiprocessing as mp

#%%
def extreme_counts(model, fixedPattern="decade", local=False):
    story = index_stats.index_stats(
        model,
        vertical_eof="ind",
        fixed_pattern=fixedPattern,
        standard="temporal_ens",
        local=local,
    )

    # temperature
    try:
        tsurf_mean = story.tsurf.mean(dim="ens").squeeze()
    except ValueError:
        tsurf_mean = story.tsurf
    tsurf_increase = tsurf_mean - tsurf_mean[0]

    # extreme counts
    ext_counts, t_surf_mean = extrc_tsurf.decadal_extrc_tsurf(
        story.eof_result.pc, tsurf_increase, ci="bootstrap"
    )

    ds = xr.Dataset({"extreme_counts": ext_counts, "tsurf": t_surf_mean})
    return ds


# %%
from multiprocessing import Pool

models = ["MPI_GE_onepct", "MPI_GE", "CanESM2", "CESM1_CAM5", "GFDL_CM3", "MK36"]


def process_model(model):
    ds = extreme_counts(model)
    ds.to_netcdf(
        f"/work/mh0033/m300883/Tel_MMLE/data/{model}/extreme_count/extre_counts_tsurf.nc"
    )


with Pool() as p:
    p.map(process_model, models)

#%%


