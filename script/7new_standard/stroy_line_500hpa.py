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
MPIGE_decade = index_stats.index_stats("MPI_GE_onepct",vertical_eof = 'ind',fixed_pattern = 'decade',standard = 'temporal_ens',season = 'MJJA',tsurf='ens_fld_year_mean')
# %%
MPIGE_decade.stat_overview(levels = np.arange(-1.8,1.9,0.3))
# %%
MPIGE_decade.extreme_count_profile(xlim = (35,110),ci = 'bootstrap')

# %%
MPIGE_decade.extrc_tsurf(ci = 'bootstrap')
#%%
MPIGE_decade.NAO_EA_hist2d()

#%%
MPIGE_decade.extrc_tsurf(ci = 'bootstrap')
#%%
MPIGE_decade.composite_analysis(reduction = 'mean_weighted')
# %%
MPIGE_decade.composite_analysis(reduction = 'mean')
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
# same for CESM1_CAM5
CESM_all = index_stats.index_stats("CESM1_CAM5",vertical_eof = 'ind',fixed_pattern = 'all',standard='temporal_ens')
CESM_all.stat_overview(levels = np.arange(-1.8,1.9,0.3))
CESM_all.extrc_tsurf(ylim = (0,60))
# %%
# same for CanESM2
Can_decade= index_stats.index_stats("CanESM2",vertical_eof = 'ind',fixed_pattern = 'decade',standard='temporal_ens')
# %%
Can_decade.stat_overview(levels = np.arange(-1.8,1.9,0.3))
# %%
Can_decade.extrc_tsurf(ylim = (0,60))
# %%
# CanESM2_all
Can_all= index_stats.index_stats("CanESM2",vertical_eof = 'ind',fixed_pattern = 'all',standard='temporal_ens')
# %%
Can_all.stat_overview(levels = np.arange(-1.8,1.9,0.3))
# %%
