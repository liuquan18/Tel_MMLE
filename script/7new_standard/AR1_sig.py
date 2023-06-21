# %%
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import statsmodels.api as sm

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


#%%
first_pc = NAO_first10.to_dataframe(name = 'NAO')['NAO'].reset_index()['NAO']
last_pc = NAO_last10.to_dataframe(name = 'NAO')['NAO'].reset_index()['NAO']
#%%
# Fit an AR1 model to the data
model_first = sm.tsa.ARIMA(first_pc, order=(1, 0, 0)).fit()
model_last = sm.tsa.ARIMA(last_pc, order=(1, 0, 0)).fit()

# %%
# Generate 9000 realizations of model_first, each 1000 long
n_realizations = 5000
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
import numpy as np
import statsmodels.api as sm

# Generate some sample data
np.random.seed(123)
nobs = 100
x = np.random.normal(size=nobs)
y = np.zeros(nobs)
for i in range(3, nobs):
    y[i] = 0.5*y[i-1] - 0.2*y[i-2] + 0.1*y[i-3] + x[i] + np.random.standard_t(3)

# Fit an AR model of order 3 with t-distributed noise
model = sm.tsa.ar_model.AutoReg(y, lags=3, trend='c', method='mle', dist='t', df=3).fit()

# Print the model summary
print(model.summary())
# %%
