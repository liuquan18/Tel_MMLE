#%%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import proplot as pplt
# %%
def read_extrc(month,model='MPI_GE_onepct'):
    dir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/extreme_count/"
    filename = f"plev_50000_decade_first_{month}_extre_counts.nc"
    ds = xr.open_dataset(dir+filename).sel(confidence = 'true')
    ds['time'] = ds.time.dt.year
    return ds.pc

def read_extrc_month(model='MPI_GE_onepct'):
    months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
    mon_idx = xr.IndexVariable('month',months)
    ds = [read_extrc(month,model) for month in months]
    ds = xr.concat(ds,dim=mon_idx)
    return ds

#%%
# MPI_GE_onepct
MPI_GE_onepct = read_extrc_month('MPI_GE_onepct')
MPI_GE_onepct_diff = MPI_GE_onepct.isel(time = -1) - MPI_GE_onepct.isel(time = 0)

# %%
# MMLEA
models = ["MPI_GE","CanESM2","CESM1_CAM5","MK36","GFDL_CM3"]
MMLEA = [read_extrc_month(model) for model in models]
#%%
# Calculate the difference between the first and last decade
MMLEA_diff = [ds.isel(time = -1) - ds.isel(time = 0) for ds in MMLEA]

#%%
# Calculate the slope of the linear regression along the 'time' dimension
MMLEA_slope = [ds.polyfit(dim='time',deg=1) for ds in MMLEA]



# %%
gs = pplt.GridSpec(nrows=2, ncols=2)
fig = pplt.figure(refwidth=2.2, span=False, share='labels')

for r, mode in enumerate(['NAO','EA']):
    for c, extr_type in enumerate(['pos','neg']):
        data = MPI_GE_onepct.sel(mode=mode,extr_type=extr_type).values
        data = np.concatenate(([data[-1]],data[:-1]),axis = 0)
        ax = fig.subplot(gs[r,c])
        lines = ax.linex(data,labels = ['first','last'])
        ax.format(title=f"{mode} {extr_type}",
                  xminorticks="null",
                  yminorlocator='null',
                  grid=False,
                  )

        ax.legend(ncols=1,loc='lr',frame=True)

        ax.set_ylim(12,-1)
        ax.set_yticks(np.arange(11,-1,-1))
        ax.set_yticklabels(["Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov"])
fig.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/monthly/month_extrc_diagram.png")

# %%
