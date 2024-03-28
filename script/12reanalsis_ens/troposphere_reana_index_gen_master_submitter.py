#%%
import numpy as np
import os
import sys
# %%
plevs = [92500, 85000, 70000, 50000, 40000, 30000, 20000]

# %%
for plev in plevs:
    os.system(f"sbatch ./troposphere_reana_index_gen_submitter.sh {plev}")
# %%
