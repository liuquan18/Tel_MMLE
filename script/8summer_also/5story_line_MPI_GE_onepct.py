#%%
import xarray as xr
import numpy as np

import src.MMLE_TEL.story_line as story_line
import src.plots.NAO_EA_hist2d as hist2d

#%%
import importlib
importlib.reload(story_line)
importlib.reload(hist2d)
# %%
MPI_GE_onepct_MJJA = story_line.story_line(
    model            =  "MPI_GE_onepct",
    vertical_eof     =  'ind',
    fixed_pattern    =  'decade',
    standard         =  'first',
    season           =   'MJJA',
    tsurf            =  'ens_fld_year_mean')
#%%

MPI_GE_onepct_MJJA.stat_overview(levels = np.arange(-1.8,1.9,0.3))
#%%
MPI_GE_onepct_MJJA.extreme_count_profile(xlim = (40,140))
#%%
MPI_GE_onepct_MJJA.extrc_tsurf()
#%%
MPI_GE_onepct_MJJA.write_doc()

#%%
MPI_GE_onepct_DJFM = story_line.story_line(
    model            =  "MPI_GE_onepct",
    vertical_eof     =  'ind',
    fixed_pattern    =  'decade',
    standard         =  'first',
    season           =   'DJFM',
    tsurf            =  'ens_fld_year_mean')

#%%
MPI_GE_onepct_DJFM.stat_overview(levels = np.arange(-1.8,1.9,0.3))
#%%
MPI_GE_onepct_DJFM.extreme_count_profile(xlim = (40,140))
#%%
MPI_GE_onepct_DJFM.extrc_tsurf()
#%%
MPI_GE_onepct_DJFM.NAO_EA_hist2d()
#%%
MPI_GE_onepct_DJFM.write_doc()


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
MPIGE_all = story_line.story_line("MPI_GE_onepct",vertical_eof = 'ind',fixed_pattern = 'all',standard='temporal_ens')
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
CESM_all = story_line.story_line("CESM1_CAM5",vertical_eof = 'ind',fixed_pattern = 'all',standard='temporal_ens')
CESM_all.stat_overview(levels = np.arange(-1.8,1.9,0.3))
CESM_all.extrc_tsurf(ylim = (0,60))
# %%
# same for CanESM2
Can_decade= story_line.story_line("CanESM2",vertical_eof = 'ind',fixed_pattern = 'decade',standard='temporal_ens')
# %%
Can_decade.stat_overview(levels = np.arange(-1.8,1.9,0.3))
# %%
Can_decade.extrc_tsurf(ylim = (0,60))
# %%
# CanESM2_all
Can_all= story_line.story_line("CanESM2",vertical_eof = 'ind',fixed_pattern = 'all',standard='temporal_ens')
# %%
Can_all.stat_overview(levels = np.arange(-1.8,1.9,0.3))
# %%
