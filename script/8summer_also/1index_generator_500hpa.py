#%%
#%%
import src.MMLE_TEL.index_generator as index_generate
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import src.compute.slurm_cluster as scluster
#%%
import importlib
importlib.reload(index_generate)
importlib.reload(scluster)

#%%


# %%
# function for generate the index
def index_gen(model,fixedPattern,plev = 50000,season = 'MJJA',standard = 'first'):
    generator = index_generate.decompose_plev(model,plev = plev,fixedPattern = fixedPattern,standard=standard,season = season)
    generator.save_result()
#%%


#%%
import concurrent.futures

with concurrent.futures.ProcessPoolExecutor() as executor:
    futures = [
        executor.submit(index_gen, 'MPI_GE_onepct', 'decade', plev=50000),
        executor.submit(index_gen, 'MPI_GE', 'decade', plev=50000),
        executor.submit(index_gen, 'CanESM2', 'decade', plev=50000),
        executor.submit(index_gen, 'CESM1_CAM5', 'decade', plev=50000),
        executor.submit(index_gen, 'MK36', 'decade', plev=50000),
        executor.submit(index_gen, 'GFDL_CM3', 'decade', plev=50000)
    ]
    for future in concurrent.futures.as_completed(futures):
        try:
            result = future.result()
        except Exception as e:
            print(f"Exception: {e}")


#%%
# MPI_GE_onepct
index_gen('MPI_GE_onepct', 'decade', plev=50000,season = 'MJJA',standard = 'first')
#%%
index_gen('MPI_GE_onepct', 'decade', plev=50000,season = 'DJFM',standard = 'first')



# %%
# CanESM2
index_gen('CanESM2', 'decade', plev=50000)

#%%
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

#%%
index_gen('MPI_GE','decade')
index_gen('MPI_GE','all')

# %%
def read_eof(model,fixedPattern):
    eof_dir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/EOF_result/plev_50000_{fixedPattern}_temporal_ens_eof_result.nc"
    eof_result = xr.open_dataset(eof_dir)
    return eof_result
#%%
# MPI_GE_onepct
MPI_GE_decade_eof = read_eof('MPI_GE_onepct','decade')

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
