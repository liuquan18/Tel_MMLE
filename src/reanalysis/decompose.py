# %%
import numpy as np

import src.MMLE_TEL.index_generator as index_generate
import src.Teleconnection.spatial_pattern as ssp

import proplot as pplt
import src.plots.statistical_overview as stat_overview
import matplotlib.pyplot as plt
# %%
from src.reanalysis.utils import read_gph_data
# %%

# %%
def decompose_period(xarr, nmode = 2,period = True):
    xarr = xarr.fillna(0) # fill nan with 0
    field = xarr.sortby("time")
    if 'ens' in xarr.dims:
        field = field.stack(com=("ens", "time"))
        dim = 'com'
    else: 
        dim = 'time'
    standard = 'eof_spatial_std' if period else 'pc_temporal_std'
    eof_result = ssp.doeof(field, standard=standard,nmode=nmode,dim = dim)
    return eof_result


# %%
def standard_period(first_eof, last_eof):
    if 'ens' in first_eof.dims:
        mean = first_eof["pc"].mean(dim=("time", "ens"))
        std = first_eof["pc"].std(dim=("time", "ens"))
    else:
        mean = first_eof["pc"].mean(dim="time")
        std = first_eof["pc"].std(dim="time")

    first_eof["pc"] = (first_eof["pc"] - mean) / std
    last_eof["pc"] = (last_eof["pc"] - mean) / std
    return first_eof, last_eof
#%%

class EOF_reanalysis:
    def __init__(self, model, group_size=40, external_forcing="quadratic_trend", period_dec = False, **kwargs):
        self.model = model
        self.group_size = group_size
        self.external_forcing = external_forcing
        self.plev = kwargs.get("plev", 50000)
        self.standard = "first"
        self.nmode = kwargs.get("nmode", 2)
        self.period_dec = period_dec
        self.variable = kwargs.get("variable", "zg")

        self.start_year = kwargs.get("start_year", "1940")
        self.end_year = kwargs.get("end_year", "2022")

        self.data = read_gph_data(model, external_forcing=external_forcing,plev = self.plev,
                                  start_year = self.start_year, end_year = self.end_year, variable = self.variable)

        if self.period_dec: # period decomposition
            self.first_data = self.data.sel(
                time=slice(self.start_year, str(int(self.start_year) + self.group_size - 1))
            )
            self.last_data = self.data.sel(
                time=slice(str(int(self.end_year) - self.group_size + 1), self.end_year)
            )

            print("********* period decomposing *********")
            self.first_eof = decompose_period(self.first_data, nmode=self.nmode)
            self.last_eof = decompose_period(self.last_data, nmode=self.nmode)
        else:
            print("********* all decomposing *********")
            self.eof = decompose_period(self.data, nmode=self.nmode, period=False)
            self.first_eof= self.eof.sel(time = slice(self.start_year, str(int(self.start_year) + self.group_size - 1)))
            self.last_eof = self.eof.sel(time = slice(str(int(self.end_year) - self.group_size + 1), self.end_year))

        print("********* standardizing *********")
        self.first_eof_std, self.last_eof_std = standard_period(
            self.first_eof, self.last_eof
        )

    def save_eof(self):
        print("********* saving *********")
        self.first_eof.to_netcdf(
            f"/work/mh0033/m300883/Tel_MMLE/data/{self.model}/EOF_result/first_{self.group_size}_eof.nc"
        )
        self.last_eof.to_netcdf(
            f"/work/mh0033/m300883/Tel_MMLE/data/{self.model}/EOF_result/last_{self.group_size}_eof.nc"
        )

        self.first_eof_std.to_netcdf(
            f"/work/mh0033/m300883/Tel_MMLE/data/{self.model}/EOF_result/first_{self.group_size}_eof_std.nc"
        )
        self.last_eof_std.to_netcdf(
            f"/work/mh0033/m300883/Tel_MMLE/data/{self.model}/EOF_result/last_{self.group_size}_eof_std.nc"
        )
        if not self.period_dec:
            self.eof.to_netcdf(
                f"/work/mh0033/m300883/Tel_MMLE/data/{self.model}/EOF_result/all_{self.group_size}_eof.nc"
            )

    def plot_spatial_index(self,levels = np.arange(-2, 2.1, 0.4), save=False):
        
        fig2 = pplt.figure(figsize=(180 / 25.4, 90 / 25.4), sharex=False, sharey=False)
        fig2.format(
            abc=True,
            abcloc="ul",
            abcstyle="a",
        )
        gs = pplt.GridSpec(
            ncols=2,
            nrows=1,
        )
        ax1 = fig2.add_subplot(gs[0], proj="ortho", proj_kw={"lon_0": -20, "lat_0": 60})
        ax2 = fig2.add_subplot(gs[1])


        if self.period_dec:
            ax1, fmap, lmap = stat_overview.spatial_pattern_plot(
                ax1,
                self.first_eof_std.eof.sel(mode="NAO").squeeze(),
                self.first_eof_std.fra.sel(mode="NAO").squeeze(),
                self.last_eof_std.eof.sel(mode="NAO", ).squeeze(),
                self.last_eof_std.fra.sel(mode="NAO").squeeze(),
                levels=levels,
            )
        else:
            ax1, fmap, lmap = stat_overview.spatial_pattern_plot(
                ax1,
                self.first_eof_std.eof.sel(mode="NAO").squeeze(),
                self.first_eof_std.fra.sel(mode="NAO").squeeze(),
                levels=levels,
            )
        ax2, hist = stat_overview.index_distribution_plot(
            ax2,
            self.first_eof_std.pc.sel(mode="NAO"),
            self.last_eof_std.pc.sel(mode="NAO"),
        )

        if save:
            plt.savefig(
                f"/work/mh0033/m300883/Tel_MMLE/docs/source/plots/supplyment/{self.model}_{self.group_size}_monthly_spatial_index.png"
            )
    def plot_ens_std(self,region = None, save = False):
        if region is not None:
            data = self.data.sel(lon = slice(region[0],region[1]), lat = slice(region[2],region[3]))
        weights = np.cos(np.deg2rad(self.data.lat))
        weights.name = "weights"
        data_weighted = self.data.weighted(weights)
        glm = data_weighted.mean(dim = ('lon','lat'))
        glm_ens_std = glm.std(dim = 'ens')
        glm_ens_std.plot.line(x = 'time')