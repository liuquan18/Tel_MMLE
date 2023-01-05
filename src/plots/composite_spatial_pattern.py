import proplot as pplt
import numpy as np
import cartopy.crs as ccrs


def composite_spatial_pattern(
    first, last, levels=np.arange(-2.0, 2.1, 0.4), hlayers=100000
):
    """
    composite map of first10 and last10 years, contourf and contour
    respectively.
    """
    fig, axes = pplt.subplots(
        space=0,
        refwidth="25em",
        wspace=3,
        hspace=3,
        proj="ortho",
        proj_kw=({"lon_0": -20, "lat_0": 60}),
        nrows=2,
        ncols=2,
    )
    axes.format(
        latlines=20,
        lonlines=30,
        coast=True,
        coastlinewidth=0.5,
        coastcolor="gray7",
        toplabels=("NAO", "EA"),
        leftlabels=("pos", "neg"),
        suptitle="Change in spatial patterns for extremes of NAO and EA",
    )

    modes = ["NAO", "EA"]  # different cols
    extr_types = ["pos", "neg"]
    for i, extr_type in enumerate(extr_types):
        for j, mode in enumerate(modes):  # one row

            first_single = first.sel(mode=mode, extr_type=extr_type, hlayers=hlayers)
            last_single = last.sel(mode=mode, extr_type=extr_type, hlayers=hlayers)

            first_m = axes[i, j].contourf(
                first_single,
                x="lon",
                y="lat",
                levels=levels,
                extend="both",
                transform=ccrs.PlateCarree(),
                cmap="RdBu_r",
            )
            axes[i, j].contour(
                last_single,
                x="lon",
                y="lat",
                color="gray8",
                nozero=True,
                labels=True,
                levels=np.delete(levels, int((len(levels) - 1) / 2)),
                labels_kw={"weight": "bold"},
            )

    fig.colorbar(first_m, loc="r", pad=3, title="gph/std")
