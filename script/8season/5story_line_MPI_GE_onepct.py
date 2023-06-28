#%%
import xarray as xr
import numpy as np

import src.MMLE_TEL.story_line as story_line
import src.plots.NAO_EA_hist2d as hist2d

#%%
import importlib
importlib.reload(story_line)
importlib.reload(hist2d)

#************************************
# %%
MPI_GE_onepct_JJAS_first = story_line.story_line(
    model            =  "MPI_GE_onepct",
    vertical_eof     =  'ind',
    fixed_pattern    =  'decade',
    standard         =  'first',
    season           =   'JJAS',
    tsurf            =  'ens_fld_year_mean')
#%%

MPI_GE_onepct_JJAS_first.stat_overview(levels = np.arange(-1.8,1.9,0.3))
#%%
MPI_GE_onepct_JJAS_first.extreme_count_profile(xlim = (40,140))
#%%
MPI_GE_onepct_JJAS_first.extrc_tsurf()
#%%
MPI_GE_onepct_JJAS_first.write_doc()

#************************************
#%%
MPI_GE_onepct_DJFM_first = story_line.story_line(
    model            =  "MPI_GE_onepct",
    vertical_eof     =  'ind',
    fixed_pattern    =  'decade',
    standard         =  'first',
    season           =   'DJFM',
    tsurf            =  'ens_fld_year_mean')

#%%
MPI_GE_onepct_DJFM_first.stat_overview(levels = np.arange(-1.8,1.9,0.3))
#%%
MPI_GE_onepct_DJFM_first.extreme_count_profile(xlim = (40,140))
#%%
MPI_GE_onepct_DJFM_first.extrc_tsurf()
#%%
MPI_GE_onepct_DJFM_first.NAO_EA_hist2d()
#%%
MPI_GE_onepct_DJFM_first.write_doc()

#************************************
# same, but for standard = 'temporal_ens'
#%%
MPI_GE_onepct_MJJA_all = story_line.story_line(
    model            =  "MPI_GE_onepct",
    vertical_eof     =  'ind',
    fixed_pattern    =  'decade',
    standard         =  'temporal_ens',
    season           =   'MJJA',
    tsurf            =  'ens_fld_year_mean')

#%%
MPI_GE_onepct_MJJA_all.extrc_tsurf()
#%%
MPI_GE_onepct_MJJA_all.extreme_count_profile(xlim = (40,140))

#%%
#************************************

MPI_GE_onepct_DJFM_all = story_line.story_line(
    model            =  "MPI_GE_onepct",
    vertical_eof     =  'ind',
    fixed_pattern    =  'decade',
    standard         =  'temporal_ens',
    season           =   'DJFM',
    tsurf            =  'ens_fld_year_mean')

#%%
MPI_GE_onepct_DJFM_all.extrc_tsurf()
#%%
MPI_GE_onepct_DJFM_all.extreme_count_profile(xlim = (40,140))
#%%
MPI_GE_onepct_DJFM_all.write_doc()

#%%
#************************************
# for season MAM with standard = 'first'
MPI_GE_onepct_MAM_first = story_line.story_line(
    model            =  "MPI_GE_onepct",
    vertical_eof     =  'ind',
    fixed_pattern    =  'decade',
    standard         =  'first',
    season           =   'MAM',
    tsurf            =  'ens_fld_year_mean')
#%%
MPI_GE_onepct_MAM_first.stat_overview(levels = np.arange(-1.8,1.9,0.3))
#%%
MPI_GE_onepct_MAM_first.extreme_count_profile(xlim = (40,140))
#%%
MPI_GE_onepct_MAM_first.extrc_tsurf()
# %%
MPI_GE_onepct_MAM_first.write_doc()
# %%

#************************************
# check JJAS but for fixed_pattern = 'all', standard = 'first', season = 'JJAS' 
MPI_GE_onepct_JJAS_all = story_line.story_line(
    model            =  "MPI_GE_onepct",
    vertical_eof     =  'ind',
    fixed_pattern    =  'all',
    standard         =  'first',
    season           =   'JJAS',
    tsurf            =  'ens_fld_year_mean')

# %%
MPI_GE_onepct_JJAS_all.stat_overview(levels = np.arange(-1.8,1.9,0.3))
# %%
MPI_GE_onepct_JJAS_all.extrc_tsurf()
# %%
