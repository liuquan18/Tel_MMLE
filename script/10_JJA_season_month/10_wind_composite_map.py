# %%
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import proplot as pplt
import src.plots.utils as utils
import matplotlib as mpl

# %%
def wind_map_single(u, v, ax,levels=np.arange(5, 21, 4)):
    # caluclate the wind speed (m/s)
    wind_speed = np.sqrt(u**2 + v**2)

    # plot the wind vector as quiver, and the wind speed as color
    contourf = ax.contourf(
        wind_speed,
        transform=ccrs.PlateCarree(), cmap="YlOrRd", 
        levels = levels,
        extend = 'both',
    )
    ax.quiver(
        u.lon[::4],
        u.lat[::4],
        u.values[::4, ::4],
        v.values[::4, ::4],
        transform=ccrs.PlateCarree(),
    )
    return ax, contourf


# %%
def wind_composite_map(U, V, plev=50000,levels = np.arange(5, 21, 4),sig = False):

    U= utils.erase_white_line(U)
    V = utils.erase_white_line(V)

    fig = pplt.figure(figsize=(150 / 25.4, 100 / 25.4),sharex=False,sharey=False)
    fig.format(
        abcloc="ul",
        abc="a",
    )
    axes = fig.subplots(
        ncols=3,
        nrows=2,
        proj="ortho",
        proj_kw=({"lon_0": -20, "lat_0": 60})
    )

    axes.format(
        latlines=20,
        lonlines=30,
        color="grey7",
        coast=True,
        coastlinewidth=0.3,
        coastcolor="charcoal",
        toplabels=["first10", "last10", "last10 - first10"],
        leftlabels=["pos", "neg"],
        suptitle=f"Change in full wind field over {plev/100:.0f} hPa",
        # set the fontsize of labels to 25
    )

    for i, extr_type in enumerate(["pos", "neg"]): # rows as extr_type
        for j, period in enumerate(["first", "last", "diff"]): # columns as period
            u = U.sel(period=period, extr_type=extr_type, mode="NAO", plev=plev)
            v = V.sel(period=period, extr_type=extr_type, mode="NAO", plev=plev)
            ax, contourf = wind_map_single(u, v, axes[i, j],levels=levels)

        # signifiance
        if sig:
            hu = axes[i,2].contourf(
                U.sel(period="diff_sig", extr_type=extr_type, mode="NAO", plev=plev),
                levels=[-0.5, 0.5, 1.5],
                colors=["none", "none"],
                hatches=["", "////"],
            )

            for c in hu.collections:
                c.set_linewidth(0.2)
                c.set_edgecolor('grey6')

            hv = axes[i,2].contourf(
                V.sel(period="diff_sig", extr_type=extr_type, mode="NAO", plev=plev),
                levels=[-0.5, 0.5, 1.5],
                colors=["none", "none"],
                hatches=["", "\\\\" ],
            )
            for c in hv.collections:
                c.set_linewidth(0.2)
                c.set_edgecolor('grey6')
    # add colorbar
    fig.colorbar(contourf, loc="r", label="Wind speed (m/s)",
                 orientation = 'vertical',shrink=0.8,width = 0.2)
    return fig
    # plt.savefig(
    #     f"/work/mh0033/m300883/Tel_MMLE/docs/source/plots/Story_line_nature_climate_change/wind_composite_map_{plev/100:.0f}hPa.png",
    # )
                 

# %%
U = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30/composite/plev_50000_decade_mpi_first_JJA_JJA_first_last_u_composite_mean.nc"
)
V = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30/composite/plev_50000_decade_mpi_first_JJA_JJA_first_last_v_composite_mean.nc"
)

#%%
wind_composite_map(U.u,V.v,plev=50000,levels = np.arange(5, 21, 4),sig = False)
plt.savefig(
    f"/work/mh0033/m300883/Tel_MMLE/docs/source/plots/Story_line_nature_climate_change/wind_composite_map_500hPa.png",
)

#%%
wind_composite_map(U.u,V.v,plev=20000,levels = np.arange(5, 21, 4),sig = False)
plt.savefig(
    f"/work/mh0033/m300883/Tel_MMLE/docs/source/plots/Story_line_nature_climate_change/wind_composite_map_200hPa.png",
)

# %%
Ur = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30/composite/plev_50000_decade_mpi_first_JJA_JJA_first_last_u_composite_mean_remove_ensmean.nc"
)
Vr = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30/composite/plev_50000_decade_mpi_first_JJA_JJA_first_last_v_composite_mean_remove_ensmean.nc"
)

#%%
wind_composite_map(Ur.u,Vr.v,plev=50000,levels = np.arange(0, 2, 0.2),sig = False)
plt.savefig(
    f"/work/mh0033/m300883/Tel_MMLE/docs/source/plots/Story_line_nature_climate_change/wind_composite_map_500hPa_remove.png",
)

#%%
fig = wind_composite_map(Ur.u,Vr.v,plev=20000,levels = np.arange(0, 2.1, 0.2),sig = False)
# change fig title
fig.suptitle(f"Change in internal wind field over 200 hPa")
plt.savefig(
    f"/work/mh0033/m300883/Tel_MMLE/docs/source/plots/Story_line_nature_climate_change/wind_composite_map_200hPa_remove.png",
)
                 
# %%
