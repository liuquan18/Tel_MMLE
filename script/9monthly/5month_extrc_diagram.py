#%%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import proplot as pplt
# %%
def read_extrc(month):
    dir = "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct/extreme_count/"
    filename = f"plev_50000_decade_first_{month}_extre_counts.nc"
    ds = xr.open_dataset(dir+filename)
    ds = ds.isel(time = [0,-1]).sel(confidence = 'true') # first and last decade
    ds['time'] = ['first','last']
    return ds.pc

# %%
months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
mon_idx = xr.IndexVariable('month',months)
# %%
ds_months = [read_extrc(month) for month in months]
ds_months = xr.concat(ds_months,dim=mon_idx)
# %%
gs = pplt.GridSpec(nrows=2, ncols=2)
fig = pplt.figure(refwidth=2.2, span=False, share='labels')

for r, mode in enumerate(['NAO','EA']):
    for c, extr_type in enumerate(['pos','neg']):
        data = ds_months.sel(mode=mode,extr_type=extr_type).values
        ax = fig.subplot(gs[r,c])
        lines = ax.linex(data,labels = ['first','last'])
        ax.format(title=f"{mode} {extr_type}",
                  xminorticks="null",
                  yminorlocator='null',
                  grid=False,
                  )

        if r == 0:
            ax.legend(ncols=1,loc='lr',frame=True)
        if r == 1:
            ax.legend(ncols=1,loc='ur',frame=True)

        ax.set_ylim(12,-1)
        ax.set_yticks(np.arange(11,-1,-1))
        ax.set_yticklabels(months)
fig.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/monthly/month_extrc_diagram.png")

# %%
