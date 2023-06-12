#%%
import os
from pathlib import Path
# %%
odir = Path("/work/mh0033/m300883/Tel_MMLE/data/")
models = ["MPI_GE_onepct", "MPI_GE", "CanESM2", "CESM1_CAM5", "GFDL_CM3", "MK36"]
for model in models:
    ts_p_dir = odir / model / "ts_processed" / "tropical_arctic_gradient.nc"
    # check if the file exists
    if os.path.exists(ts_p_dir):
        continue
    else:
        print(f"{model} doesn't exists")
# %%
