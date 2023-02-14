#%%
import src.Teleconnection.tools as tools
import src.MMLE_TEL.spatial_pattern_change as sp_change
import src.Teleconnection.spatial_pattern as ssp
import xarray as xr

# %%
data = sp_change.read_gph_data('/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/zg_processed/')
if data.hlayers.size > 1:
    data = tools.standardize(data, dim="time")

#%%
data_500 = data.sel(hlayers = 50000)

# %%
periods = [slice('1851', '1860', None),
 slice('2045', '2054', None),
 slice('2089', '2099', None)]

 
#%%
EOFs = []
FRAs = []
period_index = xr.IndexVariable(dims="period", data=periods)

# %%
for period in periods:
    print(period)
    data_p = data_500.sel(time=period)
    data_p = data_p.stack(com=("time", "ens"))
    EOF,_, FRA = ssp.doeof(data_p, nmode=2, dim="com", standard=True)
    EOFs.append(EOF)
    FRAs.append(FRA)

# %%
EOFs = xr.concat(EOFs, dim=period_index)
EOFs["period"] = ["0K", "2K", "4K"]

FRAs = xr.concat(FRAs, dim=period_index)
FRAs["period"] = ['0K', '2K', '4K']

# %%
EOFs.to_netcdf("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/EOF_result/ind_first_eofs_allPeriods.nc")
FRAs.to_netcdf("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/EOF_result/ind_first_fras_allPeriods.nc")

# %%
sp_change.spatial_pattern_maps(EOFs, FRAs)
# %%
