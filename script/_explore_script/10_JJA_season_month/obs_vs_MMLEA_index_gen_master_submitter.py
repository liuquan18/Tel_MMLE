#%%
import os
import sys

#%%
for mindex in range(6):
    os.system(f"sbatch ./obs_vs_MMLEA_index_gen_submitter.sh {mindex}") 
# %%
