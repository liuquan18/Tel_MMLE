import os
import sys
for mindex in range(6):
    os.system(f"sbatch ./osb_vs_MMLEA_index_gen_submitter.sh {mindex}") 