#%%
#%%
import src.MMLE_TEL.index_generator as index_generate
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
# %%
def standard(generator):
    pc = generator.eof_result['pc']
    std_pc = (pc - pc.mean(dim = ('time','ens')))/pc.std(dim = ('time','ens'))
    generator.std_eof_result['pc'] = std_pc
    return generator

# %%
# function for generate the index
def index_gen(model,fixedPattern):
    generator = index_generate.decompose_plev(model,plev = 50000,fixedPattern = fixedPattern,standard='temporal_ens')
    decade = standard(generator)
    decade.save_result()
# %%
# CanESM2
index_gen('CanESM2','decade')
index_gen('CanESM2','all')
# %%
# CESM1_CAM5
index_gen('CESM1_CAM5','decade')
index_gen('CESM1_CAM5','all')
# %%
# MK3.6
index_gen('MK36','decade')
index_gen('MK36','all')
#%%
# GFDL_CM3
index_gen('GFDL_CM3','decade')
index_gen('GFDL_CM3','all')


# %%
def read_eof(model,fixedPattern):
    eof_dir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/EOF_result/plev_50000_{fixedPattern}_temporal_ens_eof_result.nc"
    eof_result = xr.open_dataset(eof_dir)
    return eof_result
# %%
Can_decade_eof = read_eof('CanESM2','decade')
Can_all_eof = read_eof('CanESM2','all')
# %%
# CESM1_CAM5
CESM_decade_eof = read_eof('CESM1_CAM5','decade')
CESM_all_eof = read_eof('CESM1_CAM5','all')
# %%
# MK3.6
MK36_decade_eof = read_eof('MK36','decade')
MK36_all_eof = read_eof('MK36','all')
# %%
# GFDL_CM3
GFDL_decade_eof = read_eof('GFDL_CM3','decade')
GFDL_all_eof = read_eof('GFDL_CM3','all')
# %%
# plot the GFDL_CM3_decade EOF with the projection of projection=ccrs.Orthographic(-20,60)
def quick_plot(GFDL_decade_eof):
    GFDL_decade_eof.eof.sel(decade="1920-12-31",method = 'nearest').squeeze().plot.contourf(col = 'mode',
    subplot_kws=dict(projection=ccrs.Orthographic(-20,60)),
    levels=np.arange(-1.5,1.6,0.3),
    transform=ccrs.PlateCarree(),
)

#%%
quick_plot(GFDL_decade_eof)
# %%
quick_plot(Can_decade_eof)
# %%
quick_plot(CESM_decade_eof)
# %%
quick_plot(MK36_decade_eof)
# %%
