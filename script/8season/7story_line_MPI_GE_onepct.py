#%%
import xarray as xr
import numpy as np

import src.MMLE_TEL.story_line as story_line
import src.plots.NAO_EA_hist2d as hist2d
import src.composite.field_composite as composite


#%%
import importlib
importlib.reload(story_line)
importlib.reload(hist2d)
importlib.reload(composite)

#%%
def model_first(model,season):
    return story_line.story_line(
        model            =  model,
        vertical_eof     =  'ind',
        fixed_pattern    =  'decade',
        standard         =  'first',
        season           =   season,
        tsurf            =  'ens_fld_year_mean')

################
#%%
MPI_GE_onepct_JJAS = model_first('MPI_GE_onepct','JJAS')

#%%
MPI_GE_onepct_JJAS.stat_overview()
#%%
MPI_GE_onepct_JJAS.composite_analysis(tfield='same',levels_NAO = np.arange(-1,1.1,0.2),levels_EA = np.arange(-1,1.1,0.2))
#%%



#%%
MPI_GE_onepct_JJAS.extreme_count_profile(xlim = (40,140))
#%%
###############



#%%
def extrc_tsurf_season(season):
    model_first('MPI_GE_onepct',season).extrc_tsurf()
    model_first('MPI_GE',season).extrc_tsurf()
    model_first('CanESM2',season).extrc_tsurf()
    model_first('CESM1_CAM5',season).extrc_tsurf()
    model_first('GFDL_CM3',season).extrc_tsurf()
    model_first('MK36',season).extrc_tsurf()

#%%
extrc_tsurf_season('DJFM')

#%%
extrc_tsurf_season('MAM')



#%%
model_first('MPI_GE_onepct','JJAS').extrc_tsurf()
#%%
model_first('MPI_GE','JJAS').extrc_tsurf()
#%%
model_first('CanESM2','JJAS').extrc_tsurf()
#%%
model_first('CESM1_CAM5','JJAS').extrc_tsurf()
#%%
model_first('GFDL_CM3','JJAS').extrc_tsurf()
#%%
model_first('MK36','JJAS').extrc_tsurf()

#%%


#%%
MPI_GE_onepct_DJFM = model_first('MPI_GE_onepct','DJFM')
#%%
MPI_GE_onepct_DJFM.composite_analysis(tfield='same')
#%%
MPI_GE_onepct_DJFM.write_doc()
#%%
MPI_GE_onepct_MAM = model_first('MPI_GE_onepct','MAM')
#%%
MPI_GE_onepct_MAM.composite_analysis(tfield='same')
#%%
MPI_GE_onepct_MAM.composite_analysis(tfield='JJAS',levels_NAO=np.arange(-1,1.1,0.2),levels_EA=np.arange(-1,1.1,0.2))
#%%
MPI_GE_onepct_MAM.plot_all()
#%%
MPI_GE_onepct_MAM.extreme_count_profile(xlim = (40,140))
#%%
MPI_GE_onepct_MAM.write_doc()
#%%
MPI_GE_onepct_JJAS = model_first('MPI_GE_onepct','JJAS')

#%%
MPI_GE_onepct_JJAS.plot_all()
#%%
MPI_GE_onepct_JJAS.extreme_count_profile(xlim = (40,140))
#%%
MPI_GE_onepct_JJAS.composite_analysis(tfield='same',levels_NAO = np.arange(-1,1.1,0.2),levels_EA = np.arange(-1,1.1,0.2))

#%%
MPI_GE_onepct_JJAS.write_doc()

#%%

MPI_GE_onepct_SON = model_first('MPI_GE_onepct','SON')
MPI_GE_onepct_SON.plot_all()




#%%
MPI_GE_onepct_MAM = model_first('MPI_GE_onepct','MAM')
MPI_GE_onepct_MAM.extrc_tsurf()
#%%
MPI_GE_onepct_JJA = model_first('MPI_GE_onepct','JJA')
MPI_GE_onepct_JJA.extrc_tsurf()


#%%
CanESM2_MAM = model_first('CanESM2','MAM')
CanESM2_MAM.extrc_tsurf()

#%%
CESM1_CAM5_MAM = model_first('CESM1_CAM5','MAM')
CESM1_CAM5_MAM.extrc_tsurf()
#%%
# MK36
MK36_MAM = model_first('MK36','MAM')
MK36_MAM.extrc_tsurf()


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
