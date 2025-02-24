#%%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
# %%
def read_eof_decade(model, fixed_pattern="decade_mpi"):
    """read eofs that is decomposed by decade"""
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/EOF_result/"
    filename = f"plev_50000_{fixed_pattern}_first_JJA_eof_result.nc"
    ds = xr.open_dataset(odir + filename)
    ds = ds.sel(mode="NAO")
    return ds

# %%
eof_result = read_eof_decade("MPI_GE")
# %%
first_pc = eof_result.pc.sel(time = slice('1850','1859'))
# %%
last_pc = eof_result.pc.sel(time = slice('2090','2099'))
# %%
# 90th quantile
first_90th = first_pc.quantile(0.9) # 1.30341589

first_10th = first_pc.quantile(0.1) # -1.28656809
# %%
# Calculate the fraction of values <= 1.5
fraction_le_1_5 = (first_pc <= 1.5).mean()

# Convert fraction to percentage
percentage_le_1_5 = fraction_le_1_5 * 100 # 93%
# %%
# -1.5
fraction_le_minus_1_5 = (first_pc <= -1.5).mean()
percentage_le_minus_1_5 = fraction_le_minus_1_5 * 100 # 7%
# %%
