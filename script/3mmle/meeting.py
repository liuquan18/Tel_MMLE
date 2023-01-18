#%%
import xarray as xr
import numpy as np
import pandas as pd

import src.MMLE_TEL.index_generator as index_generate
import src.MMLE_TEL.quick_plot as quick_plot

# config
v_eof = 'ind' # vertical_eof
fpattern = 'first' # fixed pattern
# %%
cesm_qp = quick_plot.period_index("CESM1_CAM5",v_eof,fpattern,'temp')
cesm_qp.spatial_pattern_change()
# %%

## MPI-GE_onepct
mpige_onepct_qp = quick_plot.period_index("MPI_GE_onepct",v_eof,fpattern, 'temp')
# mpige_onepct_qp.plot_all()
mpige_onepct_qp.extrc_tsurf_scatter()
mpige_onepct_qp.create_doc()
# %%
## MPI
mpige_qp = quick_plot.period_index("MPI_GE",v_eof,fpattern,'temp')


# mpige_qp.plot_all()
#%%
mpige_qp.extrc_tsurf_scatter()

#%%
mpige_qp.create_doc()
# %%
## CANESM
canesm_qp = quick_plot.period_index("CanESM2",v_eof,fpattern,'temp')
canesm_qp.extrc_tsurf_scatter()
canesm_qp.create_doc()

# %%
# generate index
## CANESM2
canesm = index_generate.decompose_fixedPattern("CanESM2",v_eof,fpattern)
canesm.save_result()