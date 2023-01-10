#%%
import xarray as xr
import numpy as np
import pandas as pd

import src.MMLE_TEL.index_generator as index_generate
import src.MMLE_TEL.quick_plot as quick_plot

#%%
# config
v_eof = 'ind' # vertical_eof
fpattern = 'first' # fixed pattern

# %% 
# CANESM2
canesm = index_generate.decompose_fixedPattern("CanESM2",v_eof,fpattern)
canesm.save_result()

# %%
canesm_qp = quick_plot.period_index("CanESM2",v_eof,fpattern,'temp')
canesm_qp.plot_all()
canesm_qp.create_doc()

# %%
# MPI-GE
mpige = index_generate.decompose_fixedPattern("MPI_GE",v_eof,fpattern)
mpige.save_result()


# %%
# MPI-GE_onepct
# mpige_onepct = index_generate.decompose_fixedPattern("MPI_GE_onepct",v_eof,fpattern)
mpige_onepct_qp = quick_plot.period_index("MPI_GE_onepct",v_eof,fpattern, 'temp')
mpige_onepct_qp.plot_all()
mpige_onepct_qp.create_doc()
# %%
