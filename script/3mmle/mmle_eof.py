#%%
import xarray as xr
import numpy as np
import pandas as pd

import src.MMLE_TEL.index_generator as index_generate
import src.MMLE_TEL.quick_plot as quick_plot

# config
v_eof = 'ind' # vertical_eof
fpattern = 'decade' # fixed pattern

# %% 
# generate index
## CANESM2
canesm = index_generate.decompose_fixedPattern("CanESM2",v_eof,fpattern)

#%%
canesm.save_result()

#%%
## MPI-GE
mpige = index_generate.decompose_fixedPattern("MPI_GE",v_eof,fpattern)
mpige.save_result()

#%%
## MPI-GE_onepct
mpige_onepct = index_generate.decompose_fixedPattern("MPI_GE_onepct",v_eof,fpattern)


#%%
# NCAR mode
cesm = index_generate.decompose_fixedPattern("CESM1_CAM5",v_eof,fpattern)
cesm.save_result()


#%%
# quick plot
## CANESM
canesm_qp = quick_plot.period_index("CanESM2",v_eof,fpattern,'temp')
canesm_qp.plot_all()
canesm_qp.create_doc()

#%%
## MPI
mpige_qp = quick_plot.period_index("MPI_GE",v_eof,fpattern,'temp')


# mpige_qp.plot_all()
#%%
mpige_qp.spatial_pattern_change()
mpige_qp.extrc_tsurf_scatter()

#%%
mpige_qp.create_doc()


#%%
## MPI-GE_onepct
mpige_onepct_qp = quick_plot.period_index("MPI_GE_onepct",v_eof,fpattern, 'CO2')

#%%
# mpige_onepct_qp.plot_all()
mpige_onepct_qp.spatial_pattern_change()
mpige_onepct_qp.extrc_tsurf_scatter()
mpige_onepct_qp.create_doc()

#%%
## NCAR
cesm_qp = quick_plot.period_index("CESM1_CAM5",v_eof,fpattern,'temp')
cesm_qp.plot_all()
cesm_qp.create_doc()
# %%
