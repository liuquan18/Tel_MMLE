# %%
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import src.plots.statistical_overview as stat_overview
import numpy as np
import proplot as pplt

# %%
ERA_first = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/ERA5_allens/EOF_result/first_40_eof_std.nc"
)
EAR_last = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/ERA5_allens/EOF_result/last_40_eof_std.nc"
)
# %%
MPI_EOF = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct/EOF_result/plev_50000_decade_mpi_first_JJA_eof_result.nc"
)


# %%
def split_first_last(eof_result):
    times = eof_result.time
    years = np.unique(times.dt.year)
    first_years = years[:10]
    last_years = years[-10:]

    eof_first = eof_result.isel(decade=0).sel(
        time=eof_result["time.year"].isin(first_years)
    )
    eof_last = eof_result.isel(decade=-1).sel(
        time=eof_result["time.year"].isin(last_years)
    )
    return eof_first, eof_last


# %%
MPI_first, MPI_last = split_first_last(MPI_EOF)
# %%

fig2 = pplt.figure(figsize=(180 / 25.4, 180 / 25.4), sharex=False, sharey=False)
fig2.format(abc=True, abcloc="ul", abcstyle="a", lefttitle=["ERA5", "MPI-GE"])
gs = pplt.GridSpec(
    ncols=2,
    nrows=2,
)
ax1 = fig2.add_subplot(gs[0], proj="ortho", proj_kw={"lon_0": -20, "lat_0": 60})
ax2 = fig2.add_subplot(gs[1])
ax3 = fig2.add_subplot(gs[2], proj="ortho", proj_kw={"lon_0": -20, "lat_0": 60})
ax4 = fig2.add_subplot(gs[3])

ax1, fmap, lmap = stat_overview.spatial_pattern_plot(
    ax1,
    ERA_first.eof.sel(mode="NAO", decade="1940").squeeze(),
    ERA_first.fra.sel(mode="NAO", decade="1940").squeeze(),
    EAR_last.eof.sel(mode="NAO", decade="1982").squeeze(),
    EAR_last.fra.sel(mode="NAO", decade="1982").squeeze(),
    levels=np.arange(-2, 2.1, 0.4),
)
ax2, hist = stat_overview.index_distribution_plot(
    ax2,
    ERA_first.pc.sel(mode="NAO"),
    EAR_last.pc.sel(mode="NAO"),
)

ax3, fmap, lmap = stat_overview.spatial_pattern_plot(
    ax3,
    MPI_first.eof.sel(mode="NAO").squeeze(),
    MPI_first.fra.sel(mode="NAO").squeeze(),
    MPI_last.eof.sel(mode="NAO").squeeze(),
    MPI_last.fra.sel(mode="NAO").squeeze(),
    levels=np.arange(-2, 2.1, 0.4),
)

ax4, hist = stat_overview.index_distribution_plot(
    ax4,
    MPI_first.pc.sel(mode="NAO"),
    MPI_last.pc.sel(mode="NAO"),
)


# %%
