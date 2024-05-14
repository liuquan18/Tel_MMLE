#%% import xarray
import numpy as np
import pandas as pd
import xarray as xr
import seaborn as sns
import matplotlib.pyplot as plt

import src.plots.statistical_overview as stat_overview
import src.extreme.extreme_ci as extreme

#%%
plot2dir = '/work/mh0033/m300883/Tel_MMLE/docs/source/plots/box_based/'
doc2dir = "/work/mh0033/m300883/Tel_MMLE/docs/source/"
# %%
# read the box index
NAO = xr.open_dataset('/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct/BOX_result/NAO.nc')
EA = xr.open_dataset('/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct/BOX_result/EA.nc')

NAO = NAO.pc
EA = EA.pc
# %%
# standardize the index
NAO_std = (NAO - NAO.mean(dim = ('ens','time'))) / NAO.std(dim = ('ens','time'))
EA_std = (EA - EA.mean(dim = ('ens','time'))) / EA.std(dim = ('ens','time'))
# %%
# select the first 10 years and the last 10 years and change to df, sel the plev
def select_first10_last10(data):
    first10 = data.sel(time = slice('1851-01-01','1860-12-31'))
    last10 = data.sel(time = slice('1990-01-01','1999-12-31'))
    return first10, last10

def to_df(first10, last10,plev = 50000):
    first10 = first10.sel(plev = plev).to_dataframe().reset_index()
    first10['compare'] = ['first10']*len(first10)
    last10 = last10.sel(plev = plev).to_dataframe().reset_index()
    last10['compare'] = ['last10']*len(last10)
    both = pd.concat([first10, last10])
    both = both[['pc','compare']]
    return both

#%%
# select the first 10 years and the last 10 years
NAO_std_first10, NAO_std_last10 = select_first10_last10(NAO_std)
EA_std_first10, EA_std_last10 = select_first10_last10(EA_std)

# %%
# the distribution of the index
NAO_df = to_df(NAO_std_first10, NAO_std_last10,plev = 50000)
EA_df = to_df(EA_std_first10, EA_std_last10,plev = 50000)
sns.histplot(data = NAO_df, 
             x = 'pc', 
             hue = 'compare',
             multiple='dodge',
             shrink=.6,
             bins = np.arange(-4,4.1,0.5),
             )
plt.savefig(plot2dir + 'NAO_hist.png', dpi = 300)

#%%
sns.histplot(data = EA_df, 
             x = 'pc', 
             hue = 'compare',
             multiple='dodge',
             shrink=.6,
             bins = np.arange(-4,4.1,0.5),
             )
plt.savefig(plot2dir + 'EA_hist.png', dpi = 300)
# %%
# profile of the extreme count

## seperate pc 
mode_index = xr.IndexVariable(dims="mode", data=['NAO','EA'])
first_pc = xr.concat([NAO_std_first10, EA_std_first10], dim=mode_index)
last_pc = xr.concat([NAO_std_last10, EA_std_last10], dim=mode_index)

## count the extreme
first_count = extreme.extreme_count_xr(first_pc)
last_count = extreme.extreme_count_xr(last_pc)


# %%
extreme_profile = extreme.extreme_count_profile(
            first_count, last_count, colored=True,
            xlim = (20,100)
        )
plt.savefig(plot2dir + 'extreme_profile.png', dpi = 300)
# %%
# create md file for the plot
with open(doc2dir + 'box_based.md', 'w') as f:
    f.write('# box_based method result\n')
    f.write('## NAO distribution\n')
    f.write('![NAO distribution](plots/box_based/NAO_hist.png)\n')
    f.write('## EA distribution\n')
    f.write('![EA distribution](plots/box_based/EA_hist.png)\n')
    f.write('## Extreme profile\n')
    f.write('![Extreme profile](plots/box_based/extreme_profile.png)\n')

# %%
