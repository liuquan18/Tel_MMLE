#%%
import numpy as np
import xarray as xr
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

#%%
def read_index(model):
    file_path = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/EOF_result/plev_50000_all_first_JJA_projected_index.nc"
    index = xr.open_dataset(file_path).pc.sel(mode = 'NAO')
    return index

#%%
def cal_slope(index):
    index_notime = index.assign_coords(time = np.arange(len(index['time'])))
    # calculate the linear slope of mpi_index along time dimension
    res = index_notime.polyfit(dim = 'time', deg = 1)
    slope = res.polyfit_coefficients.sel(degree = 1)
    return slope

#%%
def convert2pd(model):
    index = read_index(model)
    slope = cal_slope(index)
    # every decade
    slope = slope * 10
    return slope.to_dataframe().reset_index()[['polyfit_coefficients']]


#%%

mpi_slope = convert2pd('MPI_GE')
canesm2_slope = convert2pd('CanESM2')
cesm1_cam5_slope = convert2pd('CESM1_CAM5')
gfdl_cm3_slope = convert2pd('GFDL_CM3')
mk36_slope = convert2pd('MK36')

#%%
# combine all dataframes into one
combined_df = pd.concat([mpi_slope, canesm2_slope, cesm1_cam5_slope, gfdl_cm3_slope, mk36_slope], axis=1)

#%%
combined_df.columns = ['MPI_GE', 'CanESM2', 'CESM1_CAM5', 'GFDL_CM3', 'MK36']
#%%

# make box plot
plt.figure(figsize=(100 / 25.4, 100 / 25.4))
sns.boxplot(data=combined_df,palette=[ "C1", "tab:purple", "tab:blue", "tab:green", "yellow"])
plt.xticks(range(5), ['MPI_GE(100)', 'CanESM2(50)', 'CESM1_CAM5(40)', 'GFDL_CM3(30)', 'MK36(20)'])
plt.ylabel('Slope / 10 yr')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/paper_supplymentary/direct_response.pdf")

# %%
