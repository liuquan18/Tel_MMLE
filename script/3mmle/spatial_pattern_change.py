#%%
import src.Teleconnection.tools as tools
import src.MMLE_TEL.spatial_pattern_change as sp_change
import src.Teleconnection.spatial_pattern as ssp
import xarray as xr
import src.warming_stage.warming_stage as warming_stage

#%%
model = "MPI_GE"
v_eof = "ind"
fpattern = "first"
odir = "/work/mh0033/m300883/Tel_MMLE/data/" + model
#%%
pc = xr.open_dataset(odir + "/EOF_result/ind_first_pc.nc")
fldmean_tsurf = warming_stage.read_tsurf_fldmean(odir + "/ts_processed/tsurf_mean.nc")
_, periods = warming_stage.split_period(pc, "temp", fldmean_tsurf)


# %%
data = sp_change.read_gph_data(
    "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/zg_processed/"
)
if data.hlayers.size > 1:
    data = tools.standardize(data, dim="time")

#%%
data_500 = data.sel(hlayers=50000)


#%%
EOFs = []
FRAs = []
period_index = xr.IndexVariable(dims="period", data=periods)

# %%
for period in periods:
    print(period)
    data_p = data_500.sel(time=period)
    data_p = data_p.stack(com=("time", "ens"))
    EOF, _, FRA = ssp.doeof(data_p, nmode=2, dim="com", standard=True)
    EOFs.append(EOF)
    FRAs.append(FRA)

# %%
EOFs = xr.concat(EOFs, dim=period_index)
EOFs["period"] = ["0K", "2K", "4K"]

FRAs = xr.concat(FRAs, dim=period_index)
FRAs["period"] = ["0K", "2K", "4K"]

# %%
EOFs.to_netcdf(
    odir+"/EOF_result/ind_first_eofs_allPeriods.nc"
)
FRAs.to_netcdf(
    odir+"/EOF_result/ind_first_fras_allPeriods.nc"
)

# %%
sp_change.spatial_pattern_maps(EOFs, FRAs)
# %%

EOFs = xr.open_dataset(
    odir+"/EOF_result/ind_first_eofs_allPeriods.nc"
).eof
FRAs = xr.open_dataset(
    odir+"/EOF_result/ind_first_fras_allPeriods.nc"
).exp_var
# %%
import matplotlib.pyplot as plt
plt.savefig("")