#%%
import xarray as xr
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
import cartopy.crs as ccrs
import seaborn as sns

from matplotlib.ticker import MaxNLocator
import matplotlib.patches as mpatches


# %%
def read_eof_rean(model, group_size=40):
    odir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/"
    first_eof_path = odir + "EOF_result/first_plev500_eof.nc"
    last_eof_path = odir + "EOF_result/last_plev500_eof.nc"

    first_eof = xr.open_dataset(first_eof_path)
    last_eof = xr.open_dataset(last_eof_path)
    return first_eof, last_eof


# read extreme counts
def read_extrc_rean(model, group_size=40):
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/extreme_count/"
    first_extc_path = odir + "first_plev50000_extc.nc"
    last_extc_path = odir + "last_plev50000_extc.nc"

    first_extc = xr.open_dataset(first_extc_path)
    last_extc = xr.open_dataset(last_extc_path)
    return first_extc, last_extc

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

def y_yerr(extrc):
    pos_true = extrc.sel(mode="NAO", confidence="true", extr_type="pos").pc.values
    pos_high = extrc.sel(mode="NAO", confidence="high", extr_type="pos").pc.values
    pos_low = extrc.sel(mode="NAO", confidence="low", extr_type="pos").pc.values
    pos_err = [pos_true - pos_low, pos_high - pos_true]

    neg_true = extrc.sel(mode="NAO", confidence="true", extr_type="neg").pc.values
    neg_high = extrc.sel(mode="NAO", confidence="high", extr_type="neg").pc.values
    neg_low = extrc.sel(mode="NAO", confidence="low", extr_type="neg").pc.values
    neg_err = [neg_true - neg_low, neg_high - neg_true]

    return pos_true, pos_err, neg_true, neg_err

def format_period_year(period, tick_number):
    if period == 0.2:
        # formater = f"first \n 1950-1979"
        formater = "first40"
    elif period == 0.8:
        # formater = f"last \n 1986-2015"
        formater = "last40"
    else:
        formater = None
        print(period)
    return formater


def reananlysis_bar(first_extrc, first_err, last_extrc, last_err, ax, 
                    x = [0.2,0.8],width = 0.4,facecolor = 'none',edgecolor = 'black',linewidth = 1,errcolor = 'black'):

    ax.bar(
        x =x,
        height = [first_extrc,last_extrc],
        width = width,
        edgecolor = edgecolor,
        facecolor = facecolor,
        linewidth = 1,
        align = 'center',
        zorder = 9,
    )
    ax.errorbar(
        x = x,
        y = [first_extrc,last_extrc],
        yerr = [first_err,last_err],
        color = errcolor,
        linewidth = 2,
        fmt='none',
        zorder = 10,
    )
    return ax
# %%
# SMILEs
fixed_pattern = "decade_mpi"
# %%
# 20CR all ens
CR20_first_extc, CR20_last_extc = read_extrc_rean("CR20_allens")

CR20_first_extc = CR20_first_extc / (4 * 80)
CR20_last_extc = CR20_last_extc / (4 * 80)

# %%
CR20_ens_first_eof, CR20_ens_last_eof = read_eof_rean("CR20")

CR20_ens_first_pc = CR20_ens_first_eof.pc.sel(mode="NAO") - CR20_ens_first_eof.pc.sel(mode = 'NAO').mean()
CR20_ens_last_pc = CR20_ens_last_eof.pc.sel(mode="NAO") - CR20_ens_last_eof.pc.sel(mode = 'NAO').mean()

# also read ensemble mean of 20CR
CR20_ens_first_extc, CR20_ens_last_extc = read_extrc_rean("CR20")
CR20_ens_first_extc = CR20_ens_first_extc / 4
CR20_ens_last_extc = CR20_ens_last_extc / 4


#%%
CR20_composite_ts_first = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/composite/first_composite_mean_ts_40.nc")
CR20_composite_ts_last = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/composite/last_composite_mean_ts_40.nc")
CR20_composite_ts_diff = CR20_composite_ts_last - CR20_composite_ts_first
CR20_composite_ts = xr.concat([CR20_composite_ts_first, CR20_composite_ts_last, CR20_composite_ts_diff], dim="period")
CR20_composite_ts = CR20_composite_ts.ts
CR20_composite_ts['period'] = ["first", "last", "diff"]

#%%
CR20_composite_psl_first = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/composite/first_composite_mean_psl_40.nc")
CR20_composite_psl_last = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/composite/last_composite_mean_psl_40.nc")
CR20_composite_psl_diff = CR20_composite_psl_last - CR20_composite_psl_first
#%%
# from Pa to hPa
CR20_composite_psl_first = CR20_composite_psl_first.ts / 100
CR20_composite_psl_last = CR20_composite_psl_last.ts / 100
CR20_composite_psl_diff = CR20_composite_psl_diff.ts / 100


#%%
pc_first_df = CR20_ens_first_eof.pc.sel(mode="NAO").to_dataframe().reset_index()
pc_last_df = CR20_ens_last_eof.pc.sel(mode="NAO").to_dataframe().reset_index()

#%%
pc_first_df["period"] = "first"
pc_last_df["period"] = "last"

pc_dfs = pd.concat([pc_first_df[['pc', 'period']], pc_last_df[['pc', 'period']]], ignore_index=True)


#%%
fig = plt.figure(figsize=(12, 8))

# Only map axes (2 rows x 3 columns)
ax_pos_first = fig.add_subplot(231, projection=ccrs.Orthographic(central_longitude=-20, central_latitude=60))
ax_pos_last = fig.add_subplot(232, projection=ccrs.Orthographic(central_longitude=-20, central_latitude=60))
ax_pos_diff = fig.add_subplot(233, projection=ccrs.Orthographic(central_longitude=-20, central_latitude=60))

ax_neg_first = fig.add_subplot(234, projection=ccrs.Orthographic(central_longitude=-20, central_latitude=60))
ax_neg_last = fig.add_subplot(235, projection=ccrs.Orthographic(central_longitude=-20, central_latitude=60))
ax_neg_diff = fig.add_subplot(236, projection=ccrs.Orthographic(central_longitude=-20, central_latitude=60))


# positive
# ts as color
CR20_composite_ts_first.ts.sel(extr_type = 'pos', mode = 'NAO').plot.contourf(
    ax=ax_pos_first,
    transform=ccrs.PlateCarree(),
    cmap="RdBu_r",
    levels =np.arange(-1.5, 1.6, 0.3),
    add_colorbar=False,
)

CR20_composite_ts_last.ts.sel(extr_type = 'pos', mode = 'NAO').plot.contourf(
    ax=ax_pos_last,
    transform=ccrs.PlateCarree(),
    cmap="RdBu_r",
    levels =np.arange(-1.5, 1.6, 0.3),
    extend = 'both',
    add_colorbar=False, 
)
CR20_composite_ts_diff.ts.sel(extr_type = 'pos', mode = 'NAO').plot.contourf(
    ax=ax_pos_diff,
    transform=ccrs.PlateCarree(),
    cmap="RdBu_r",
    levels =np.arange(-1.5, 1.6, 0.3),
    extend = 'both',
    add_colorbar=False,
)

# psl as lines
CR20_composite_psl_first.sel(extr_type = 'pos', mode = 'NAO').plot.contour(
    ax=ax_pos_first,
    transform=ccrs.PlateCarree(),
    colors="black",
    linewidths =1,
    add_colorbar=False,
    levels = [l for l in np.arange(-5, 5.1, 1) if l != 0],
)

CR20_composite_psl_last.sel(extr_type = 'pos', mode = 'NAO').plot.contour(
    ax=ax_pos_last,
    transform=ccrs.PlateCarree(),
    colors="black",
    linewidths =1,
    add_colorbar=False,
    levels = [l for l in np.arange(-5, 5.1, 1) if l != 0],
)
CR20_composite_psl_diff.sel(extr_type = 'pos', mode = 'NAO').plot.contour(
    ax=ax_pos_diff,
    transform=ccrs.PlateCarree(),
    colors="black",
    linewidths =1,
    add_colorbar=False, 
    levels = [l for l in np.arange(-5, 5.1, 1) if l != 0],
)
# negative

map_temp = CR20_composite_ts_first.ts.sel(extr_type='neg', mode='NAO').plot.contourf(
    ax=ax_neg_first,
    transform=ccrs.PlateCarree(),
    cmap="RdBu_r",
    levels=np.arange(-1.5, 1.6, 0.3),
    extend='both',
    add_colorbar=False,
)

CR20_composite_ts_last.ts.sel(extr_type='neg', mode='NAO').plot.contourf(
    ax=ax_neg_last,
    transform=ccrs.PlateCarree(),
    cmap="RdBu_r",
    levels=np.arange(-1.5, 1.6, 0.3),
    extend='both',
    add_colorbar=False,
)

CR20_composite_ts_diff.ts.sel(extr_type='neg', mode='NAO').plot.contourf(
    ax=ax_neg_diff,
    transform=ccrs.PlateCarree(),
    cmap="RdBu_r",
    levels=np.arange(-1.5, 1.6, 0.3),
    extend='both',
    add_colorbar=False,
)

CR20_composite_psl_first.sel(extr_type='neg', mode='NAO').plot.contour(
    ax=ax_neg_first,
    transform=ccrs.PlateCarree(),
    colors="black",
    linewidths =1,
    add_colorbar=False,
    levels=[l for l in np.arange(-5, 5.1, 1) if l != 0],
)

CR20_composite_psl_last.sel(extr_type='neg', mode='NAO').plot.contour(
    ax=ax_neg_last,
    transform=ccrs.PlateCarree(),
    colors="black",
    linewidths =1,
    add_colorbar=False,
    levels=[l for l in np.arange(-5, 5.1, 1) if l != 0],
)

CR20_composite_psl_diff.sel(extr_type='neg', mode='NAO').plot.contour(
    ax=ax_neg_diff,
    transform=ccrs.PlateCarree(),
    colors="black",
    linewidths =1,
    add_colorbar=False,
    levels=[l for l in np.arange(-5, 5.1, 1) if l != 0],
)


for ax in [ax_pos_first, ax_pos_last, ax_pos_diff, ax_neg_first, ax_neg_last, ax_neg_diff]:
    ax.set_global()
    ax.coastlines(linewidth=0.5)
    ax.gridlines(linestyle="--")
    ax.set_title("")




# Add colorbar to the right of the plots
cax = fig.add_axes([1.0, 0.3, 0.015, 0.5])  # [left, bottom, width, height]

fig.colorbar(
    map_temp,
    cax=cax,
    orientation="vertical",
    aspect=30,
    label="near surface temperature anomaly (K)",
)

plt.tight_layout()

# # save the figure
fig.savefig("/work/mh0033/m300883/Tel_MMLE/docs/defense/20CR_allens_nao_extreme_maps_cb.png",
            bbox_inches="tight",
            dpi=500,
            )

# %%
