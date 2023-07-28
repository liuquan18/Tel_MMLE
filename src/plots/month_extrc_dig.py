#%%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import proplot as pplt
import matplotlib.patches as mpatches
# %%
def bar_hatch(xlim=(-0.4, 0.6)):

    gs = pplt.GridSpec(nrows=1, ncols=2)
    fig = pplt.figure(refwidth=2.2, span=False, share="labels")
    # the right order of the models
    models = ["MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]
    models_legend = [
        "MPI-GE (100)",
        "CanESM2 (50)",
        "CESM1-CAM5 (40)",
        "MK3.6 (30)",
        "GFDL-CM3 (20)",
    ]
    colors_model = ["C1", "tab:purple", "tab:blue", "tab:green", "C4"]
    model_color = dict(zip(models, colors_model))
    months = [
        "Nov",
        "Oct",
        "Sep",
        "Aug",
        "Jul",
        "Jun",
        "May",
        "Apr",
        "Mar",
        "Feb",
        "Jan",
        "Dec",
    ]

    for r, mode in enumerate(["NAO"]):
        for c, extr_type in enumerate(["pos", "neg"]):
            ax = fig.subplot(gs[r, c])
            ax.format(
                title=f"{mode} {extr_type}",
                xminorticks="null",
                yminorlocator="null",
                grid=False,
                xlim=xlim,
                ylim=(12, -1),
                yticks=np.arange(11, -1, -1),
                yticklabels=months,
            )
            for i, month in enumerate(months):
                data_month = MMLEA_slope.sel(
                    mode=mode, extr_type=extr_type, month=month
                )
                data_cum = data_month.cumsum(axis=0).values
                hatch_cum = abs(data_month.values).cumsum(axis=0)
                for j, model in enumerate(models):
                    data_model = data_month.sel(model=model).values
                    bar_width = data_model
                    starts = data_cum[j] - bar_width
                    hatch = None
                    zorder = 0
                    # Check if the current model's data has a different sign than the previous model's data
                    if (j > 0 and data_model * data_month[j - 1] < 0) or (
                        j > 1
                        and any([abs(data_cum[j]) < abs(data_cum[i]) for i in range(j)])
                    ):
                        # If either of these conditions is true, set the hatch pattern for the current bar to "/////"
                        hatch = "/////"
                    bar = ax.barh(
                        i,
                        bar_width,
                        fc=colors_model[j],
                        edgecolor="blue9",
                        label=models_legend[j],
                        left=starts,
                        hatch=hatch,
                        width=1.6,
                        alpha=1,
                        zorder=zorder,
                    )

    patchs = [
        mpatches.Patch(facecolor=colors_model[i], label=models[i], edgecolor="blue9")
        for i in range(len(models))
    ]
    patchs.append(
        mpatches.Patch(
            facecolor="none", label="overlap", hatch="/////", edgecolor="blue9"
        )
    )
    fig.legend(patchs, ncols=3, loc="b", frame=False)

    fig.savefig(
        "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/monthly/month_extrc_diagram_patch.png"
    )


# %%
def bar_side(xlim=(-0.3, 0.3)):
    gs = pplt.GridSpec(nrows=1, ncols=2)
    fig = pplt.figure(refheight=4.2, refwidth=2.2, span=False, share="labels")
    # the right order of the models
    models = ["MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]
    models_legend = [
        "MPI-GE (100)",
        "CanESM2 (50)",
        "CESM1-CAM5 (40)",
        "MK3.6 (30)",
        "GFDL-CM3 (20)",
    ]
    colors_model = ["C1", "tab:purple", "tab:blue", "tab:green", "C4"]
    model_color = dict(zip(models, colors_model))
    months = [
        "Nov",
        "Oct",
        "Sep",
        "Aug",
        "Jul",
        "Jun",
        "May",
        "Apr",
        "Mar",
        "Feb",
        "Jan",
        "Dec",
    ]

    for r, mode in enumerate(["NAO"]):
        for c, extr_type in enumerate(["pos", "neg"]):
            data = MMLEA_slope.sel(mode=mode, extr_type=extr_type)
            y_pos = np.arange(0, 12, 1)
            bar_width = 0.15

            ax = fig.subplot(gs[r, c])
            ax.format(
                title=f"{mode} {extr_type}",
                xminorticks="null",
                yminorlocator="null",
                grid=False,
                xlim=xlim,
                ylim=(12, -1),
                yticks=y_pos + 0.37,
                yticklabels=months,
            )

            for i, model in enumerate(MMLEA_slope.model.values):
                ax.barh(
                    y_pos + i * bar_width,
                    data.sel(model=model).values,
                    bar_width,
                    fc=colors_model[i],
                    edgecolor="blue9",
                    label=models_legend[i],
                )

    patchs = [
        mpatches.Patch(facecolor=colors_model[i], label=models[i], edgecolor="blue9")
        for i in range(len(models))
    ]
    patchs.append(
        mpatches.Patch(
            facecolor="none", label="overlap", hatch="/////", edgecolor="blue9"
        )
    )
    fig.legend(patchs, ncols=3, loc="b", frame=False)
    fig.savefig(
        "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/monthly/month_extrc_diagram_patch_side.png"
    )
