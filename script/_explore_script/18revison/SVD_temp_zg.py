# %%
import logging.config
import xarray as xr
import numpy as np
from xmca.array import MCA
from xmca.xarray import MCA as xMCA
import glob
import logging
from src.Teleconnection.tools import sqrtcoslat
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
# info logging
logging.basicConfig(level=logging.INFO)
# %%
# read data
def read_data_month(var, month, name = None, plev = None):
    if name is None:
        name = var
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/{var}_{month}/"
    files = glob.glob(odir + "*.nc")
    # sort
    files.sort()
    # read
    data_month = xr.open_mfdataset(files, combine='nested', concat_dim='ens', parallel=True)
    data_month = data_month[name]
    if plev is not None:
        try:
            data_month = data_month.sel(plev = plev)
        except KeyError:
            logging.warning(f"Variable {name} does not have plev dimension")
            pass
    return data_month

def read_data(var, name = None, plev = 50000, remove_ens_mean = False, temporal_shift = 0):
    if name is None:
        name = var
    data = []
    for month in ['May', 'Jun','Jul','Aug']: # include may to shift
        data_month = read_data_month(var, month, name, plev)
        data.append(data_month)

    data = xr.concat(data, dim='time')    
    data = data.sortby('time')

    if temporal_shift != 0:
        data = data.shift(time = temporal_shift)

    data = data.sel(time = data.time.dt.month.isin([6,7,8])) # only keep three months

    # remove ensemble mean
    if remove_ens_mean:
        data = data - data.mean(dim='ens')

    first = data.sel(time = slice('1850','1859'))
    last = data.sel(time = slice('2090', '2099'))
    if var == 'mrso':
        last = last.isel(time = slice(0, 30)) # 2100 is labelled as 2099

    first = first.stack(com=('ens','time')).transpose('com','lat','lon')
    last = last.stack(com=('ens','time')).transpose('com','lat','lon')

    return first.compute(), last.compute()

#%%
def eofs_to_xarray(data, eof, pc, fra):

    modes = np.arange(pc.shape[-1])

    reduce_dim = data.dims[0] # 'com' or 'time'
    eof_cnt = data[0]
    eof_cnt = eof_cnt.drop_vars(reduce_dim) # drop the dim 'com' or 'time'
    eof_cnt = eof_cnt.expand_dims(dim = {'mode':modes},axis = 0) # add the dim 'mode' and 'decade'

    eof_cnt = eof_cnt.transpose('lat','lon', 'mode')
    eofx = eof_cnt.copy(data=eof)

    # pc to xarray
    pcx = xr.DataArray(
        pc, dims=[reduce_dim, "mode"], coords={reduce_dim: data[reduce_dim], "mode":modes}
    )
    pcx = pcx.unstack()
    frax = xr.DataArray(fra, dims=["mode",], coords={"mode": modes})
    return eofx,pcx,frax


#%%
zg_first, zg_last = read_data('zg', 'var156')

# %%
ts_first, ts_last = read_data('ts', 'tsurf', None, temporal_shift = 1)
# %%
mrso_first, mrso_last = read_data('mrso', 'var140', None)


#%%
def do_mca(west, east, nmodes = 5):
    west = (west - west.mean(dim='com')) / west.std(dim='com')
    east = (east - east.mean(dim='com')) / east.std(dim='com')

    west = west.fillna(0)
    east = east.fillna(0)

    west_w = sqrtcoslat(west)
    east_w = sqrtcoslat(east)

    west = west * west_w
    east = east * east_w

 
    mca = MCA(west.values, east.values)
    mca.solve()

    exp_var = mca.explained_variance(nmodes)
    pcs = mca.pcs(nmodes)
    eofs = mca.eofs(nmodes)

    return exp_var, pcs, eofs
#%%
def mca_eof(west, east):
    exp_var, pcs, eofs = do_mca(west, east, nmodes=5)

    lef_pcs, right_pcs = [pcs['left'], pcs['right']]

    left_eofs, right_eofs = [eofs['left'], eofs['right']]

    left_eofs, lef_pcs, lef_fra = eofs_to_xarray(west, left_eofs, lef_pcs, exp_var)
    right_eofs, right_pcs, right_fra = eofs_to_xarray(east, right_eofs, right_pcs, exp_var)

    # de weight
    left_weight = sqrtcoslat(left_eofs)
    right_weight = sqrtcoslat(right_eofs)

    left_eofs = left_eofs / left_weight
    right_eofs = right_eofs / right_weight

    return left_eofs, lef_pcs, lef_fra, right_eofs, right_pcs, right_fra

#%%
zg_zt_eofs_first, zg_zt_pcs_first, zg_zt_fra_first, ts_zt_eofs_first, ts_zt_pcs_first, ts_zt_fra_first = mca_eof(zg_first, ts_first)
zg_zt_eofs_last, zg_zt_pcs_last, zg_zt_fra_last, ts_zt_eofs_last, ts_zt_pcs_last, ts_zt_fra_last = mca_eof(zg_last, ts_last)

# %%
mrso_mt_eofs_first, mrso_mt_pcs_first, mrso_mt_fra_first, ts_mt_eofs_first, ts_mt_pcs_first, ts_mt_fra_first = mca_eof(mrso_first, ts_first)
mrso_mt_eofs_last, mrso_mt_pcs_last, mrso_mt_fra_last, ts_mt_eofs_last, ts_mt_pcs_last, ts_mt_fra_last = mca_eof(mrso_last, ts_last)

#%%
zg_first_std = zg_first.std(dim='com')
ts_first_std = ts_first.std(dim='com')
mrso_first_std = mrso_first.std(dim='com')

zg_last_std = zg_last.std(dim='com')
ts_last_std = ts_last.std(dim='com')
mrso_last_std = mrso_last.std(dim='com')

# back to original scale
zg_zt_eofs_first = zg_zt_eofs_first * zg_first_std 
ts_zt_eofs_first = ts_zt_eofs_first * ts_first_std
mrso_mt_eofs_first = mrso_mt_eofs_first * mrso_first_std 

zg_zt_eofs_last = zg_zt_eofs_last * zg_last_std 
ts_zt_eofs_last = ts_zt_eofs_last * ts_last_std
mrso_mt_eofs_last = mrso_mt_eofs_last * mrso_last_std

#%%
ts_levels = np.arange(-0.1, 0.11, 0.02)
zg_levels = np.arange(-1.5, 1.51, 0.3)*2
mr_levels = np.arange(-0.15, 0.151, 0.03)/50
#%%
fig, axes = plt.subplots(2,3, figsize = (12,8), subplot_kw={'projection': ccrs.Orthographic(-20, 60)})
zg_zt_eofs_first.sel(mode=3).plot(ax = axes[0,0], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = zg_levels, extend = 'both', cbar_kwargs={'label':'Z500','shrink':0.8})
(zg_zt_eofs_last.sel(mode=3)*-1).plot(ax = axes[0,1], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = zg_levels, extend = 'both', cbar_kwargs={'label':'Z500','shrink':0.8})
((zg_zt_eofs_last.sel(mode=3)*-1) - zg_zt_eofs_first.sel(mode=3)).plot(ax = axes[0,2], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = zg_levels, extend = 'both', cbar_kwargs={'label':'Z500','shrink':0.8})

ts_zt_eofs_first.sel(mode=3).plot(ax = axes[1,0], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = ts_levels, extend = 'both', cbar_kwargs={'label':'temperature','shrink':0.8})
(ts_zt_eofs_last.sel(mode=3)*-1).plot(ax = axes[1,1], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = ts_levels, extend = 'both', cbar_kwargs={'label':'temperature','shrink':0.8})
((ts_zt_eofs_last.sel(mode=3)*-1) - ts_zt_eofs_first.sel(mode=3)).plot(ax = axes[1,2], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = ts_levels, extend = 'both', cbar_kwargs={'label':'temperature','shrink':0.8})

for ax in axes.flat:
    ax.coastlines()
    ax.set_global()
    ax.set_title('')

# add a, b, c
axes[0,0].text(0.05, 0.9, 'a', transform=axes[0,0].transAxes, fontsize=16, fontweight='bold')
axes[0,1].text(0.05, 0.9, 'b', transform=axes[0,1].transAxes, fontsize=16, fontweight='bold')
axes[0,2].text(0.05, 0.9, 'c', transform=axes[0,2].transAxes, fontsize=16, fontweight='bold')
axes[1,0].text(0.05, 0.9, 'd', transform=axes[1,0].transAxes, fontsize=16, fontweight='bold')
axes[1,1].text(0.05, 0.9, 'e', transform=axes[1,1].transAxes, fontsize=16, fontweight='bold')
axes[1,2].text(0.05, 0.9, 'f', transform=axes[1,2].transAxes, fontsize=16, fontweight='bold')

plt.tight_layout()
# plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/paper_supplymentary/zg_ts_SVD.pdf")
plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/paper_supplymentary/zg_ts_1monshift_SVD.pdf")

#%%
fig, axes = plt.subplots(2,3, figsize = (12,8), subplot_kw={'projection': ccrs.Orthographic(-20, 60)})
(mrso_mt_eofs_first.sel(mode=3)*-1).plot(ax = axes[0,0], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = mr_levels, extend = 'both', cbar_kwargs={'label':'soil wetness','shrink':0.8})
(mrso_mt_eofs_last.sel(mode=3)*1).plot(ax = axes[0,1], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = mr_levels, extend = 'both', cbar_kwargs={'label':'soil wetness','shrink':0.8})
((mrso_mt_eofs_last.sel(mode=3)*1) - (mrso_mt_eofs_first.sel(mode=3)*-1)).plot(ax = axes[0,2], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = mr_levels, extend = 'both', cbar_kwargs={'label':'soil wetness','shrink':0.8})

(ts_mt_eofs_first.sel(mode=3)*-1).plot(ax = axes[1,0], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = ts_levels, extend = 'both', cbar_kwargs={'label':'temperature','shrink':0.8})
(ts_mt_eofs_last.sel(mode=3)*1).plot(ax = axes[1,1], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = ts_levels, extend = 'both', cbar_kwargs={'label':'temperature','shrink':0.8})
((ts_mt_eofs_last.sel(mode=3)*1) - (ts_mt_eofs_first.sel(mode=3)*-1)).plot(ax = axes[1,2], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = ts_levels, extend = 'both', cbar_kwargs={'label':'temperature','shrink':0.8})

for ax in axes.flat:
    ax.coastlines()
    ax.set_global()
    ax.set_title('')

# add a, b, c
axes[0,0].text(0.05, 0.9, 'a', transform=axes[0,0].transAxes, fontsize=16, fontweight='bold')
axes[0,1].text(0.05, 0.9, 'b', transform=axes[0,1].transAxes, fontsize=16, fontweight='bold')
axes[0,2].text(0.05, 0.9, 'c', transform=axes[0,2].transAxes, fontsize=16, fontweight='bold')
axes[1,0].text(0.05, 0.9, 'd', transform=axes[1,0].transAxes, fontsize=16, fontweight='bold')
axes[1,1].text(0.05, 0.9, 'e', transform=axes[1,1].transAxes, fontsize=16, fontweight='bold')
axes[1,2].text(0.05, 0.9, 'f', transform=axes[1,2].transAxes, fontsize=16, fontweight='bold')

plt.tight_layout()
# plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/paper_supplymentary/mrso_ts_SVD.pdf")
plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/paper_supplymentary/mrso_ts_1monshift_SVD.pdf")





# %%
fig, axes = plt.subplots(2,4, figsize = (12,12), subplot_kw={'projection': ccrs.Orthographic(-20, 60)})
zg_zt_eofs_first.sel(mode=4).plot(ax = axes[0,0], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = zg_levels, extend = 'both')
zg_zt_eofs_first.sel(mode=1).plot(ax = axes[0,1], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = zg_levels, extend = 'both')
zg_zt_eofs_first.sel(mode = 2).plot(ax = axes[0,2], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = zg_levels, extend = 'both')
zg_zt_eofs_first.sel(mode = 3).plot(ax = axes[0,3], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = zg_levels, extend = 'both')

zg_zt_eofs_last.sel(mode=4).plot(ax = axes[1,0], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = zg_levels, extend = 'both')
zg_zt_eofs_last.sel(mode=1).plot(ax = axes[1,1], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = zg_levels, extend = 'both')
zg_zt_eofs_last.sel(mode = 2).plot(ax = axes[1,2], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = zg_levels, extend = 'both')
zg_zt_eofs_last.sel(mode = 3).plot(ax = axes[1,3], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = zg_levels, extend = 'both')


for ax in axes.flat:
    ax.coastlines()
    ax.set_global()
    ax.set_title('')
plt.show()

#%%
fig, axes = plt.subplots(2,4, figsize = (12,12), subplot_kw={'projection': ccrs.Orthographic(-20, 60)})
ts_zt_eofs_first.sel(mode=0).plot(ax = axes[0,0], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = ts_levels, extend = 'both')
ts_zt_eofs_first.sel(mode=1).plot(ax = axes[0,1], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = ts_levels, extend = 'both')
ts_zt_eofs_first.sel(mode = 2).plot(ax = axes[0,2], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = ts_levels, extend = 'both')
ts_zt_eofs_first.sel(mode = 3).plot(ax = axes[0,3], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = ts_levels, extend = 'both')

ts_zt_eofs_last.sel(mode=0).plot(ax = axes[1,0], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = ts_levels, extend = 'both')
ts_zt_eofs_last.sel(mode=1).plot(ax = axes[1,1], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = ts_levels, extend = 'both')
ts_zt_eofs_last.sel(mode = 2).plot(ax = axes[1,2], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = ts_levels, extend = 'both')
ts_zt_eofs_last.sel(mode = 3).plot(ax = axes[1,3], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = ts_levels, extend = 'both')


for ax in axes.flat:
    ax.coastlines()
    ax.set_global()
    ax.set_title('')
plt.show()
# %%

mt_exp_var, mt_pcs, mt_eofs = do_mca(mrso_last, ts_last, nmodes=5)
# %%
mrso_mt_pcs, ts_mt_pcs = [mt_pcs['left'], mt_pcs['right']]
mrso_mt_eofs, ts_mt_eofs = [mt_eofs['left'], mt_eofs['right']]

# %%
mrso_mt_eofs, mrso_mt_pcs, mrso_mt_fra = eofs_to_xarray(mrso_first, mrso_mt_eofs, mrso_mt_pcs, mt_exp_var)
ts_mt_eofs, ts_mt_pcs, ts_mt_fra = eofs_to_xarray(ts_first, ts_mt_eofs, ts_mt_pcs, mt_exp_var)
# %%
fig, axes = plt.subplots(2,3, figsize = (12,10), subplot_kw={'projection': ccrs.Orthographic(-20, 60)})
mrso_mt_eofs.sel(mode=0).plot(ax = axes[0,0], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = mr_levels, extend = 'both')
mrso_mt_eofs.sel(mode=3).plot(ax = axes[0,1], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = mr_levels, extend = 'both')
mrso_mt_eofs.sel(mode = 4).plot(ax = axes[0,2], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = mr_levels, extend = 'both')

ts_mt_eofs.sel(mode=0).plot(ax = axes[1,0], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = t_levels, extend = 'both')
ts_mt_eofs.sel(mode=3).plot(ax = axes[1,1], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = t_levels, extend = 'both')
ts_mt_eofs.sel(mode=4).plot(ax = axes[1,2], transform=ccrs.PlateCarree(), cmap='coolwarm', levels = t_levels, extend = 'both')

for ax in axes.flat:
    ax.coastlines()
    ax.set_global()
    ax.set_title('')
plt.show()

# %%
