#%%
import xarray as xr
import numpy as np
import src.Teleconnection.spatial_pattern as ssp
# %%
model = 'CESM1_CAM5'


# %%
def correct_sign(model):
    dir = '/work/mh0033/m300883/Tel_MMLE/data/' + model + '/EOF_result/ind_first_'

    print("read data...")
    eofx = xr.open_dataset(dir+'old_eof.nc').eof
    pcx = xr.open_dataset(dir+'old_pc.nc').pc

    # change sign
    print("correct sign")
    coef = ssp.sign_coef(eofx)
    eofx = eofx * coef
    pcx = pcx * coef

    eofx.name = 'eof'
    pcx.name = 'pc'

    eofx.to_netcdf(dir+'eof.nc')
    pcx.to_netcdf(dir+'pc.nc')
# %%
correct_sign('CESM1_CAM5')
# %%
correct_sign('CanESM2')
# %%
