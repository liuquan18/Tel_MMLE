#%%
import os
from multiprocessing import Pool
from pathlib import Path
import numpy as np
# %%
models = ["MPI_GE_onepct", "MPI_GE", "CanESM2", "CESM1_CAM5", "GFDL_CM3", "MK36"]
odir = Path("/work/mh0033/m300883/Tel_MMLE/data/")

for model in models:
    for plev in ['50000','30000']:
        # try if the file exists
        file_dir = odir / model / "extreme_count" / f"plev_{plev}_extre_counts_tsurf.nc"
        if not os.path.exists(file_dir):
            continue

        # raname the file name to "plev_{plev}_temporal_ens_extre_counts_tsurf.nc"
        new_file_dir = odir / model / "extreme_count" / f"plev_{plev}_temporal_ens_extre_counts_tsurf.nc"
        os.rename(file_dir,new_file_dir)
# %%
for ens_size in np.arange(20,101,10):
    random_file_dir = odir / "MPI_GE_onepct_random" / "extreme_count" / f"plev_30000_extreme_counts_tsurf_{ens_size}.nc"
    random_new_file_dir = odir / "MPI_GE_onepct_random" / "extreme_count" / f"plev_30000_temporal_ens_extreme_counts_tsurf_{ens_size}.nc"
    os.rename(random_file_dir,random_new_file_dir)
# %%
for ens_size in np.arange(20,101,10):
    random_file_dir = odir / "MPI_GE_onepct_random" / "extreme_count" / f"plev_50000_extreme_counts_tsurf_{ens_size}.nc"
    random_new_file_dir = odir / "MPI_GE_onepct_random" / "extreme_count" / f"plev_50000_temporal_ens_extreme_counts_tsurf_{ens_size}.nc"
    os.rename(random_file_dir,random_new_file_dir)

# %%
