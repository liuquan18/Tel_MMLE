# %%
import xarray as xr
import pandas as pd
import seaborn as sns
# %%

def read_extrc(model, fixed_pattern="decade_mpi"):
    """read extreme counts"""
    odir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/extreme_count/"
    filename = f"plev_50000_{fixed_pattern}_first_JJA_extre_counts.nc"
    ds = xr.open_dataset(odir + filename).pc

    # divide the ensemble size of each model
    ens_sizes = {
        "MPI_GE": 100,
        "MPI_GE_onepct": 100,
        "CanESM2": 50,
        "CESM1_CAM5": 40,
        "MK36": 30,
        "GFDL_CM3": 20,
    }
    ds = ds / ens_sizes[model]
    count = ds.sel(confidence = 'true').sel(mode = 'NAO')
    return count

# %%
extreme_counts = read_extrc('MPI_GE')
# %%
NAO = extreme_counts.to_dataframe('count').reset_index()[['time','extr_type', 'count']]
# %%
uex = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/u/hist_rcp85_0002_echam6_ATM_mm.nc")
# %%
