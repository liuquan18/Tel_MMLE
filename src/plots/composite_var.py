import cartopy.crs as ccrs
import proplot as pplt
import src.plots.utils as utils


def composite_var(var, first, last, mode, unit=None):

    if var == "tsurf":
        unit = "K"

    first = utils.erase_white_line(first)
    last = utils.erase_white_line(last)

    data_all = [
        first.sel(mode=mode),
        last.sel(mode=mode),
        last.sel(mode=mode) - first.sel(mode=mode),
    ]
    extr_type = ["pos", "neg"]

    fig, axes = pplt.subplots(
        space=0,
        refwidth="25em",
        wspace=3,
        hspace=3,
        proj="ortho",
        proj_kw=({"lon_0": -20, "lat_0": 60}),
        nrows=2,
        ncols=3,
    )
    axes.format(
        latlines=20,
        lonlines=30,
        coast=True,
        coastlinewidth=0.5,
        coastcolor="gray7",
        toplabels=["first10", "last10", "last10 - first10"],
        leftlabels=("pos", "neg"),
        suptitle=f"Change in influence of extreme {mode} on {var}",
    )

    extr_types = ["pos", "neg"]
    for i, extr_type in enumerate(extr_types):
        for j, data in enumerate(data_all):  # one row

            first_m = axes[i, j].contourf(
                data.sel(extr_type=extr_type),
                x="lon",
                y="lat",
                # levels=levels,
                extend="both",
                transform=ccrs.PlateCarree(),
                cmap="RdBu_r",
            )

    fig.colorbar(first_m, loc="r", pad=3, title=f"{var}/{unit}")
