#%%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
# %%
def pre_process(ds):
    return ds.sortby("time")
def read_data(model,var = 'ts'):
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/{var}_"
    try:
        f_Jun = xr.open_mfdataset(f"{odir}Jun/*.nc", combine="by_coords")
        f_Jul = xr.open_mfdataset(f"{odir}Jul/*.nc", combine="by_coords")
        f_Aug = xr.open_mfdataset(f"{odir}Aug/*.nc", combine="by_coords")
    except ValueError:
        f_Jun = xr.open_mfdataset(f"{odir}Jun/*.nc", combine="nested", concat_dim="ens")
        f_Jul = xr.open_mfdataset(f"{odir}Jul/*.nc", combine="nested", concat_dim="ens")
        f_Aug = xr.open_mfdataset(f"{odir}Aug/*.nc", combine="nested", concat_dim="ens")
    f = xr.concat([f_Jun, f_Jul, f_Aug], dim="time")
    f = f.sortby("time")
    return f

# %%
def glmean(xarr):
    weights = np.cos(np.deg2rad(xarr.lat))
    xarr_weighted = xarr.weighted(weights)
    return xarr_weighted.mean(("lon", "lat"))

# anomaly to 1940-1979
def _anomaly(xarr):
    return xarr - xarr.sel(time=slice("1979", "2000")).mean("time")
def anomaly(xarr):
    return xarr.groupby("time.month").apply(_anomaly)

# rolling mean
def _rolling_mean(xarr, window = 8):
    rm = xarr.rolling(time = window, center = True).mean()
    return rm
def rolling_mean(xarr, window = 8):
    rm = xarr.groupby("time.month").apply(_rolling_mean, window = window)
    return rm

# quadratic fit 
def poly_fit(xarr):
    xarr_yearly = xarr.copy()
    xarr_yearly = xarr_yearly.resample(time = '1Y').mean()
    poly_coef = xarr_yearly.polyfit(dim="time", deg=2)
    poly_fit = xr.polyval(xarr.time, poly_coef)
    return poly_fit
def _poly_fit(xarr):
    coef = xarr.polyfit(dim="time", deg=2)
    poly_fit = xr.polyval(xarr.time, coef)
    return poly_fit
def quadratic_fit(xarr):
    poly_fit = xarr.groupby("time.month").apply(_poly_fit)
    return poly_fit
#%%
def plot_month_series(ERA5_mean, Allens_mean, ens_mean, ens_mean_rm, ens_mean_fit, month = 7,period = slice('1979','2022')):
    if period is not None:
        ERA5_mean = ERA5_mean.sel(time = period)
        Allens_mean = Allens_mean.sel(time = period)
        ens_mean = ens_mean.sel(time = period)
        ens_mean_rm = ens_mean_rm.sel(time = period)
        ens_mean_fit = ens_mean_fit.sel(time = period)
    fig, ax = plt.subplots(figsize=(10, 5))
    ERA5_mean.sel(time = ERA5_mean.time.dt.month == month).plot(ax=ax, label="ERA5",zorder = 3)
    ens_mean.sel(time = ens_mean.time.dt.month == month).plot(ax = ax, label = 'ensmean',zorder = 4)
    ens_mean_rm.sel(time = ens_mean_rm.time.dt.month == month).plot(ax = ax, label = 'ensmean_rm',zorder = 5)
    ens_mean_fit.sel(time = ens_mean_fit.time.dt.month == month).plot(ax = ax, label = 'ensmean_fit',zorder = 6,color = 'r')
    [Allens.plot(ax=ax,color = 'grey',zorder = 1) for Allens in Allens_mean.sel(time = Allens_mean.time.dt.month == month)]
    # ax.set_ylim(-0.5,1.25)
    ax.legend()


def plot_series(ERA5_mean, Allens_mean, ens_mean, ens_mean_rm, ens_mean_fit, period = slice('1979','2022')):
    if period is not None:
        ERA5_mean = ERA5_mean.sel(time = period)
        Allens_mean = Allens_mean.sel(time = period)
        ens_mean = ens_mean.sel(time = period)
        ens_mean_rm = ens_mean_rm.sel(time = period)
        ens_mean_fit = ens_mean_fit.sel(time = period)

    fig, ax = plt.subplots(figsize=(10, 5))
    ERA5_mean.plot(ax=ax, label="ERA5",zorder = 3)
    ens_mean.plot(ax = ax, label = 'ensmean',zorder = 4)
    ens_mean_rm.plot(ax = ax, label = 'ensmean_rm',zorder = 5)
    ens_mean_fit.plot(ax = ax, label = 'ensmean_fit',zorder = 6,color = 'r')
    [Allens.plot(ax=ax,color = 'grey',zorder = 1) for Allens in Allens_mean]
    # ax.set_ylim(-0.5,1.25)
    ax.legend()
# %%
# temperature
TERA5 = read_data("ERA5", var = 'zg')

#%%
TAllens = read_data("CR20_allens")
# %%
ERA5_mean = glmean(TERA5)
Allens_mean = glmean(TAllens)
#%%
# ERA5_mean = ERA5_mean.sel(time = slice('1979','2022'))
# Allens_mean = Allens_mean.sel(time = slice('1979','2022'))

# %%
ERA5_mean = anomaly(ERA5_mean.T2M)
Allens_mean = anomaly(Allens_mean.ts)
ens_mean = Allens_mean.mean(dim = 'ens')

ens_mean_rm = rolling_mean(ens_mean, window = 8)
# %%
ens_mean_fit = poly_fit(ens_mean_rm).polyfit_coefficients

# %%
plot_month_series(ERA5_mean, Allens_mean, ens_mean, ens_mean_rm, ens_mean_fit)
plot_series(ERA5_mean, Allens_mean, ens_mean, ens_mean_rm, ens_mean_fit)
#%%
Tres = Allens_mean - ens_mean_fit


# %%
# geopotential height
ZERA5 = read_data("ERA5",var = 'zg')
ZERA5 = ZERA5.sel(plev = 50000)
ZERA5 = ZERA5/9.81
Zallens = read_data("ERA5_allens",var = 'zg')
# %%
ZERA5_s = glmean(ZERA5.var129.sel(lat = slice(60,35), lon = slice(-30,0)))
Zallens_s = glmean(Zallens.zg.sel(lat = slice(60,35), lon = slice(-30,0)))
# %%
ZERA5_n = glmean(ZERA5.var129.sel(lat = slice(75,60), lon = slice(-70,-50)))
Zallens_n = glmean(Zallens.zg.sel(lat = slice(75,60), lon = slice(-70,-50)))
# %%
ZERA5_s = anomaly(ZERA5_s)
Zallens_s = anomaly(Zallens_s)
ZERA5_n = anomaly(ZERA5_n)
Zallens_n = anomaly(Zallens_n)
# %%
Zens_mean_s = Zallens_s.mean(dim = 'ens')
Zens_mean_n = Zallens_n.mean(dim = 'ens')
# %%
ZERA5_s_rm = rolling_mean(ZERA5_s, window = 8)
Zens_mean_s_rm = rolling_mean(Zens_mean_s, window = 8)
# %%
ZERA5_n_rm = rolling_mean(ZERA5_n, window = 8)
Zens_mean_n_rm = rolling_mean(Zens_mean_n, window = 8)

# %%
# fit quadratic
# Zens_mean_s_fit = poly_fit(Zens_mean_s_rm).polyfit_coefficients
# Zens_mean_n_fit = poly_fit(Zens_mean_n_rm).polyfit_coefficients

Zens_mean_s_fit = quadratic_fit(Zens_mean_s_rm).polyfit_coefficients
Zens_mean_n_fit = quadratic_fit(Zens_mean_n_rm).polyfit_coefficients

# %%
ZERA5_n_fit = _poly_fit(ZERA5_n_rm).polyfit_coefficients
ZERA5_s_fit = _poly_fit(ZERA5_s_rm).polyfit_coefficients
# %%
plot_month_series(ZERA5_s, Zallens_s, Zens_mean_s, Zens_mean_s_rm, Zens_mean_s_fit, month = 7)

# %%
plot_series(ZERA5_s, Zallens_s, Zens_mean_s, Zens_mean_s_rm, Zens_mean_s_fit)
# %%

#%%
def _iv(xarr):
    iv = xarr.std()
    return iv
IV_era5 = Tres.resample(time = '1Y').apply(_iv)
# %%
# %%
Zres_s = Zallens_s - Zens_mean_s_fit
# %%
IV_s = Zres_s.resample(time = '20Y').apply(_iv)

# %%
