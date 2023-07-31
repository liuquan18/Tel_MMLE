#%%
import xarray as xr
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import proplot as pplt

# %%
def read_extrc(model):
    odir = '/work/mh0033/m300883/Tel_MMLE/data/'+model+'/extreme_count/'
    filename = 'plev_50000_decade_mpi_first_JJA_extre_counts.nc'
    ds = xr.open_dataset(odir+filename).pc
    ds = ds.sel(confidence = 'true')
    return ds
# %%
models = ['MPI_GE','CanESM2','CESM1_CAM5','MK36','GFDL_CM3']
extrcs = [read_extrc(model) for model in models]
# %%
gs = pplt.GridSpec(nrows=1, ncols=2)
fig = pplt.figure(refwidth=2.2, refheight = 5.2, span=False, share="labels")
# the right order of the models
models = ["MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]
models_legend = [
    "MPI-GE (100)",
    "CanESM2 (50)",
    "CESM1-CAM5 (40)",
    "MK3.6 (30)",
    "GFDL-CM3 (20)",
]
colors_model = ["C1", "tab:purple", "tab:blue", "tab:green", "C4"]
model_color = dict(zip(models, colors_model))

lines = []
for r, mode in enumerate(['NAO']):
    for c, extr_type in enumerate(['pos','neg']):
        ax = fig.subplot(gs[r, c])
        for model in models:
            extrc = extrcs[models.index(model)]
            line = extrc.sel(mode=mode,extr_type=extr_type).plot(
                ax=ax, 
                label=model_color[model],
                x = 'time',
                color = model_color[model],
                linewidth = 2,)
            lines.append(line)
            
        ax.format(
            ylim=(20, 280),
            ylabel="Extreme counts",
            xlabel="Year",
            title=f"{mode} {extr_type}",
            suptitle="",
            titleloc="uc",
            ylocator=20,
            yminorlocator="null",
            grid=False,

        )
fig.legend(
    lines,
    labels=models_legend,
    ncols=3,
    loc="b",
)
plt.savefig(
    '/work/mh0033/m300883/Tel_MMLE/docs/source/plots/monthly/JJA_extreme_counts_line_time_all.png',
)
# %%
