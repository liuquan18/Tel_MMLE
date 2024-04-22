#%%
import numpy as np
import os
import sys
# %%
# plevs = [92500, 85000, 70000, 50000, 40000, 30000, 20000]

plevs = [1000, 975, 950, 925, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200]
model="CR20"
# %%
for plev in plevs:
    os.system(f"sbatch ./troposphere_reana_index_gen_submitter.sh {plev} {model}")
# %%
