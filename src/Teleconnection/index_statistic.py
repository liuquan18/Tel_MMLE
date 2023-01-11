"""
source code for index statistics
"""
import numpy as np
import pandas as pd
import seaborn as sns
import xarray as xr
from eofs.standard import Eof
from tqdm.notebook import tqdm, trange

import src.Teleconnection.temporal_index as sti


def df_sel(df,condition):
    """
    select dataframe rows depend on the condition of items in columns. {'column':'item'}
    """
    for item in condition.items():
        col = item[0]
        value = item[1]
        df = df[df[col]==value]

    return df

def pc_todf(xarr,name):
    """
    transform xarray to dataframe.
    **Arguments**
        *xarr* the dataarray to be transoformed
        *name* the column name for the dataframe.
    """
    df = xarr.__xarray_dataarray_variable__.to_dataframe().reset_index()
    df = df.rename(columns = {'__xarray_dataarray_variable__':name})
    return df

def merge_df(dfs,on=['hlayers','mode','ens','time']):
    """
    merge dataframe along the selected columns. 
    """

    dfm = pd.merge(dfs[0],dfs[1],on = on)

    dfm = pd.merge(dfm,dfs[2],on = on)

    old_columns = dfs[0].columns[:4]
    new_columns = [dfs[0].columns[-1],dfs[1].columns[-1],dfs[2].columns[-1]]

    dfm.columns = np.r_[old_columns,new_columns]
    return dfm

def show_allyears_ens(df,mode):
    Z200_EA = df_sel(df,{'hlayers':20000,'mode':mode}).set_index(['hlayers','mode','ens','time'])
    Z500_EA = df_sel(df,{'hlayers':50000,'mode':mode}).set_index(['hlayers','mode','ens','time'])
    Z850_EA = df_sel(df,{'hlayers':85000,'mode':mode}).set_index(['hlayers','mode','ens','time'])
    return Z200_EA,Z500_EA,Z850_EA


def pc_column(pcs,mode='NAO'):
    """
    combine all three pcs make them into three columns of dataframe.
    **Arguments**
        *pcs* list of pcs in xarray, [first,all,last]
        *mode* mode to plot.
    """
    dfs = [pc.sel(mode=mode).to_dataframe()[['pc']] for pc in pcs]
    df = dfs[0].join(dfs[1],lsuffix = '_first',rsuffix = '_all')
    df = df.join(dfs[2])
    df.columns = ['pc_first','pc_all','pc_last']

    df = df.stack()
    df = pd.DataFrame(df,columns = ['pc'])
    df.index.names = ['hlayers','ens','time','spap']
    return df

def join_xarr(lxarr, rxarr,lsuffix,rsuffix = "_all"):
    """
    join two xarry into dataframes, with the first column being the index 
    from first(last)-pattern, and the second the index from all-pattern.
    **Arguments**
        *lxarr* the index of ten period onto first or last pattern.
        *rxarr* the index of ten period onto all pattern
        *lsuffix* "_first" or "_last"
        *rsuffix* "_all"
    **Return**
        one dataframe with the first column as "pc_first" or "pc_last" (lsuffix).
        the second column as "pc_all"
    """
    ten_ten_all = lxarr.to_dataframe().join(rxarr.to_dataframe(),
    lsuffix = lsuffix, rsuffix = rsuffix)  # ten_ten, ten_all
    return ten_ten_all

def pattern_compare(all_indexes):
    """
    compare the index from first (last) pattern and all pattern. 
    **Argument**
        *all_indexes* index of projecting all the time series onto three patterns.
        *period* "first","last","dynamic", corresponding to the comparation between
            - first-first --- first-all, first-last --- first-all
            - last-first  --- last-all,  last-last  --- last-all
            - first-first --- first-all, last-last --- last-all
            here the first 'first' represent period, the second 'first' represent pattern.
    **Return**
        dataframe with the first column being the index from first pattern ,the second
        column being the index from all.
    """
    # getting the index for the two periods
    first_all, first_first,first_last = sti.period_index(all_indexes,period = 'first10') # ten-->period
    last_all, last_first, last_last = sti.period_index(all_indexes,period = 'last10')
    
    # first 10 periods
    first_first_all = join_xarr(first_first,first_all,lsuffix='_first')
    first_last_all  = join_xarr(first_last, first_all,lsuffix='_last')

    # last 10 periods
    last_first_all = join_xarr(last_first,last_all,lsuffix='_first')
    last_last_all = join_xarr(last_last, last_all,lsuffix='_last')

    return [first_first_all, first_last_all],[last_first_all,last_last_all],\
        [first_first_all,last_last_all]


#################### extreme counts #########################
def all_layer_counts(extc):
    """
    sum the counts of all hlayers together.
    """
    extc = extc.reset_index()
    extc_all = extc.groupby(['extr_type','mode','pattern','period'])[['extreme_counts']].sum()
    extc_all = extc_all.reset_index().set_index(['extr_type','mode'])
    return extc_all

def count_nonzero(xarr):
    """
    count the number of extreme cases in xarr
    **Arguments**
        *xarr* the xarr where the non-extreme points are marked as np.nan
    **Return**
        number of extreme cases, with the coordinate of 'hlayer' and 'mode'
        reserved.
    """
    return xarr.count(dim = ('time','ens'))

def extreme_count(all_indexes):
    """
    count the number of extreme cases in the period of first 10 years and 
    last 10 years.
    **Arguments**
        *all_indexes* the three index of all years from the three patterns.
                      should be order in [all_all, all_first, all_last]
    **Return**
        three for each period (six in total) extreme counts.
        *first10_ec* the extreme count of first 10 period.
        *last10_ec* the extreme count of last 10 period.
    """
    # getting the extreme indexes
    first10_extremes = sti.period_extreme(all_indexes,period = 'first10') # first 10 period
    last10_extremes = sti.period_extreme(all_indexes, period = 'last10') # last 10 period

    # count the extreme cases
    first10_ec = [count_nonzero(first10_e) for first10_e in first10_extremes]
    last10_ec = [count_nonzero(last10_e) for last10_e in last10_extremes]

    return first10_ec,last10_ec


def extr_count_df(first_extreCounts,last_extreme_Counts):
    """
    transform xarray to dataframe. two columns: first 10 and last 10 period.
    new index level pattern: first, all and last.
    **Arguments**
        *extreCounts* [fA,fF,fL,lA,lF,lL]
        *hlayer* which vertical layer to be transform.
    **Return**
        dataframe, with two columns and a new index level "pattern"
    """
    # first 10 period
    patterns = ['all','first','last']
    first_df = [first_exCount.to_dataframe(name = patterns[i])[[patterns[i]]]
    for i, first_exCount in enumerate(first_extreCounts)]
    first_df = pd.concat(first_df,axis = 1)
    first_df = pd.DataFrame(first_df.stack(),columns = ['first10'])

    # last 10 period
    last_df = [last_exCount.to_dataframe(name = patterns[i])[[patterns[i]]]
    for i, last_exCount in enumerate(last_extreme_Counts)]
    last_df = pd.concat(last_df,axis = 1)
    last_df = pd.DataFrame(last_df.stack(),columns=['last10'])

    # final df
    final_df = first_df.join(last_df)
    final_df = pd.DataFrame(final_df.stack(),columns=['extreme_counts'])
    final_df.index.names = ['extr_type','hlayers','mode','pattern','period']
    final_df = final_df.reset_index().set_index(['extr_type','mode','hlayers'])

    return final_df


def period_diff(extreme_count):
    """
    calculate the difference of the number of extreme events between the first 10 and
    the last 10 years.
    **Arguments**
        *extreme_count* the dataframe of extreme counts, gotten from function 
        'extr_count_df'
    **Return**
        *diff* diff dataframe, with four columns representing different patterns.
    """

    # first pattern data 
    fF = df_sel(extreme_count,{'period':'first10','pattern':'first'})
    lF = df_sel(extreme_count,{'period':'last10', 'pattern':'first'})

    # last pattern data
    fL = df_sel(extreme_count,{'period':'first10','pattern':'last' })
    lL = df_sel(extreme_count,{'period':'last10', 'pattern':'last'})

    # all pattern data
    fA = df_sel(extreme_count,{'period':'first10','pattern':'all'})
    lA = df_sel(extreme_count,{'period':'last10', 'pattern':'all'})

    # first pattern diff
    firstDiff = lF[['extreme_counts']]-fF[['extreme_counts']]

    # last pattern diff
    lastDiff = lL[['extreme_counts']]-fL[['extreme_counts']]

    # dynamic diff
    dynamicDiff = lL[['extreme_counts']]-fF[['extreme_counts']]

    # all diff
    allDiff = lA[['extreme_counts']]-fA[['extreme_counts']]

    # together
    diff = firstDiff.join(lastDiff,lsuffix = '_first',rsuffix = '_last')
    diff = diff.join(dynamicDiff)
    diff = diff.join(allDiff,lsuffix = '_dynamic',rsuffix = '_all')


    # column names
    columns = pd.Index(['first','last','dynamic','all'],name = 'diff')
    diff.columns = columns
    return diff


def combine_diff(extreme_counts,mode):
    """
    combine the first10, last 10 and the difference between those two into one
    dataframe.
    **Arguments**
        *extreme_count* the dataframe of extreme counts. from function 'extr_count_df'.
    **Return**
        *all* list of dataframes [first10,last10,diff]
    """
    diff = period_diff(extreme_counts) # get the data for thir panel
    diff = diff.xs(mode, level = 'mode')
    extreme_counts = extreme_counts.xs(mode,level = 'mode')

    # unstack diff
    diff = diff.unstack(0)

    # unstack extreme counts
    extreme_counts = extreme_counts.reset_index()\
    .set_index(['hlayers','pattern','extr_type','period'])\
        .unstack([1,2,3])

    # split to first10 and last 10
    first10 = extreme_counts.xs('first10',level = 'period',axis=1)['extreme_counts']
    last10 = extreme_counts.xs('last10',level = 'period',axis=1)['extreme_counts']

    all = [first10,last10,diff]
    return all

def extc_diff_first(extc):
    """
    different of rise of extreme events from the first pattern.
    **Argument**
        *extc* the increment of extrc cases.
    **Return**
        *diff* the difference between lastpattern/allpattern and firstpattern as dataframe
    """
    AF = extc['all']-extc['first']
    LF = extc['last']-extc['first']

    # dataframe
    diff = pd.DataFrame([AF,LF])
    diff = diff.T
    diff.columns = ['AF','LF']

    # put 'AF' and 'LF' as index
    diff = pd.DataFrame(diff.stack(),columns=['extr_diff'])

    # put 'pos' and 'neg' as columns
    diff = diff.unstack(0)['extr_diff']
    diff = diff.reset_index().rename(columns={'level_2':'diff_pattern'})
    diff = diff.set_index(['diff_pattern','hlayers','mode'])
    return diff

def uncenter_corr(
    a: np.array,
    b: np.array
)-> np.float32:
    """
    the uncentered correlation of two lists
    """
    return np.dot(a,b)/np.dot(b,b)


def xr_corr(
    spa_A: xr.DataArray,
    spa_B: xr.DataArray,
    dims: str,
    center: bool = True
)-> xr.DataArray:
    """
    calculate the spatial corrlation of two dataArray.
    calculated use both uncenter and centered pattern correlations to 
    assess the spatial similarity with and without accounting for the
    magnitude of the patterns respectively.
    **Arguments**
        *spa_A* the data A.
        *spa_B* the data B.
        *dims* along which to calculate correlation.
        *center* center or uncenter.
    **Return**
        *Corr* the correlation data Array.
    """
    if center:
        Corr = xr.corr(spa_A,spa_B,dim = dims)
    else:
        spa_A = spa_A.stack(com = dims)
        spa_B = spa_B.stack(com = dims)
        Corr = xr.apply_ufunc(  uncenter_corr,
                                spa_A,
                                spa_B,
                                vectorize=True,
                                dask='parallelized',
                                input_core_dims=[['com'],['com']],
        
        )
    return Corr

def spatial_corr(
    All: xr.DataArray,
    First: xr.DataArray,
    Last: xr.DataArray,
    center: bool = True
)-> pd.DataFrame:
    """
    similarity of all to first and last to first.
    calculate the spatial similarity of all to first and last to first 
    with (center=False) and without (center = True) accounting for the
    magnitude of the pattern.
    **Arguments**
        *spa_A* the data A.
        *spa_B* the data B.
        *dims* along which to calculate correlation.
        *center* center or uncenter.
    **Return**
        *Corr* the correlation dataframe.
    """
    AF = xr_corr(All,First,dims = ['lat','lon'],center=center)
    LF = xr_corr(Last,First, dims = ['lat','lon'],center=center)

    index = xr.DataArray(['AF','LF'],dims = 'diff_pattern')
    corr = xr.concat([AF,LF],dim = index)
    corr = corr.to_dataframe()[['eof']]
    return corr




