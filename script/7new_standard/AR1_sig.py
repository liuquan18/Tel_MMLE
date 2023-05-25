# %%
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
# %%
eof_result = xr.open_dataset('/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct/EOF_result/plev_50000_decade_temporal_ens_eof_result.nc')
# %%
pc = eof_result['pc']
# %%
NAO = pc.sel(mode='NAO')

# %%
NAO_first10 = NAO.isel(time=slice(0, 10))
NAO_last10 = NAO.isel(time=slice(-10, None))

# %%
NAO_first10 = NAO_first10.stack(realization=('ens', 'time'))
NAO_last10 =  NAO_last10.stack(realization=('ens', 'time'))
# %%
NAO_first_count = (NAO_first10 > 1.5).sum().item()
NAO_last_count = (NAO_last10 > 1.5).sum().item()
# %%
import numpy as np

# Define a function to calculate the count of values above a threshold
def count_above_threshold(data, threshold = 1.5):
    return (data > threshold).sum().item()

# Define a function to generate a bootstrap sample
def bootstrap_sample(data):
    return np.random.choice(data, size=len(data))

# Define a function to generate bootstrap replicates
def bootstrap_replicates(data, func, n_reps):
    replicates = np.empty(n_reps)
    for i in range(n_reps):
        sample = bootstrap_sample(data)
        replicates[i] = func(sample)
    return replicates

# Calculate the bootstrap replicates of NAO_first_count and NAO_last_count
n_reps = 10000
NAO_first_replicates = bootstrap_replicates(NAO_first10, count_above_threshold, n_reps)
NAO_last_replicates = bootstrap_replicates(NAO_last10, count_above_threshold, n_reps)

# Calculate the confidence intervals
conf_int_first = np.percentile(NAO_first_replicates, [2.5, 97.5])
conf_int_last = np.percentile(NAO_last_replicates, [2.5, 97.5])

# Print the confidence intervals
print(f"The 95% confidence interval for NAO_first_count is [{conf_int_first[0]:.2f}, {conf_int_first[1]:.2f}].")
print(f"The 95% confidence interval for NAO_last_count is [{conf_int_last[0]:.2f}, {conf_int_last[1]:.2f}].")
# %%
import statsmodels.api as sm

#%%
first_pc = NAO_first10.to_dataframe(name = 'NAO')['NAO'].reset_index()['NAO']
last_pc = NAO_last10.to_dataframe(name = 'NAO')['NAO'].reset_index()['NAO']
#%%
# Fit an AR1 model to the data
model_first = sm.tsa.ARIMA(first_pc, order=(1, 0, 0)).fit()
model_last = sm.tsa.ARIMA(last_pc, order=(1, 0, 0)).fit()

# Calculate the confidence interval using the AR1 model
conf_int_first = model_first.conf_int(alpha=0.05)
conf_int_last = model_last.conf_int(alpha=0.05)

# Print the confidence interval
print(f"The 95% confidence interval for NAO_first_count is [{conf_int_first.iloc[1, 0]:.2f}, {conf_int_first.iloc[1, 1]:.2f}].")
print(f"The 95% confidence interval for NAO_last_count is [{conf_int_last.iloc[1, 0]:.2f}, {conf_int_last.iloc[1, 1]:.2f}].")
# %%
one_real = model_first.predict(start = 0, end = 999)
# %%
# Generate 9000 realizations of model_first, each 1000 long
n_realizations = 9000
n_obs = 1000
simulations = np.empty((n_realizations, n_obs))
for i in range(n_realizations):
    simulations[i, :] = model_first.simulate(nsimulations=n_obs)

# Print the shape of the simulations array
print(f"The shape of the simulations array is {simulations.shape}.")
# %%
# count the number of values above 1.5 in simulations along the second axis
counts = (simulations > 1.5).sum(axis=1)
# %%
# calculate the 5th and 95th percentiles of counts
percentiles = np.percentile(counts, [5, 95])
# %%
