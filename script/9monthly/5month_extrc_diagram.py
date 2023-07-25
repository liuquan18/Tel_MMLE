#%%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import proplot as pplt
# %%
def read_extrc(month,model='MPI_GE_onepct'):
    dir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/extreme_count/"
    filename = f"plev_50000_decade_first_{month}_extre_counts.nc"
    ds = xr.open_dataset(dir+filename)
    ds = ds.isel(time = [0,-1]).sel(confidence = 'true') # first and last decade
    ds['time'] = ['first','last']
    return ds.pc

# %%
months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
mon_idx = xr.IndexVariable('month',months)
# %%
# MPI_GE_onepct
MPI_GE_onepct = [read_extrc(month) for month in months]
MPI_GE_onepct = xr.concat(MPI_GE_onepct,dim=mon_idx)

# %%
# MMLEA
MMLEA_diff = []
for model in ['CanESM2','CESM1_CAM5','MK36','GFDL_CM3']:
    for month in months:
        first_last = read_extrc(month,model)
        MMLEA_diff.append(first_last.sel(time='last') - first_last.sel(time='first'))

# %%
# count how many of the 5 models, the last - first is positive



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
