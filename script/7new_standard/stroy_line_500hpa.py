#%%
import xarray as xr
import numpy as np

import src.MMLE_TEL.index_stats as index_stats
import src.plots.NAO_EA_hist2d as hist2d

#%%
import importlib
importlib.reload(index_stats)
importlib.reload(hist2d)
# %%
MPIGE_decade = index_stats.index_stats("MPI_GE_onepct",vertical_eof = 'ind',fixed_pattern = 'decade',standard = 'temporal_ens')
# %%
MPIGE_decade.stat_overview(levels = np.arange(-1.8,1.9,0.3))
# %%
MPIGE_decade.extrc_tsurf()
#%%
MPIGE_decade.NAO_EA_hist2d()
#%%
MPIGE_decade.write_doc()
# %%

# same for fixedpattern = 'all'
MPIGE_all = index_stats.index_stats("MPI_GE_onepct",vertical_eof = 'ind',fixed_pattern = 'all',standard='temporal_ens')
# %%
MPIGE_all.stat_overview(levels = np.arange(-1.8,1.9,0.3))
# %%
MPIGE_all.extrc_tsurf()
#%%
MPIGE_all.NAO_EA_hist2d()
#%%
MPIGE_all.write_doc()
# %%
