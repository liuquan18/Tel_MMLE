#%%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import proplot as pplt
import matplotlib.patches as mpatches

# %%
def bar_hatch(MMLEA_slope, xlim=(-0.35, 0.65)):

    gs = pplt.GridSpec(nrows=1, ncols=2)
    fig = pplt.figure(refwidth=2.2, span=False, share="labels")
    # the right order of the models
    models = MMLEA_slope.model.values
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
        "Dec",
        "Jan",
        "Feb",
        "Mar",
        "Apr",
        "May",
        "Jun",
        "Jul",
        "Aug",
        "Sep",
        "Oct",
        "Nov",
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
                yticks=np.arange(12),
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
                    zorder = j
                    # the first block, shorter than the second block, and the sign changed,
                    # hatch the first block
                    # leave the second block as it is
                    if (
                        j == 0
                        and data_model * data_month[j + 1] < 0
                        and abs(data_model) < abs(data_month[j + 1])
                    ):
                        hatch = "/////"
                        zorder = j + 1

                    # the second or later block, shorter than the previous block, and the sign changed,
                    # hatch the current block
                    if (
                        j > 0
                        and data_model * data_month[j - 1] < 0
                        and abs(data_model) < abs(data_month[j - 1])
                    ):
                        hatch = "/////"

                    # there is one  block that is longer than any of the previous blocks
                    # hatch the current block
                    if any([abs(data_cum[j]) < abs(data_cum[i]) for i in range(j)]):
                        hatch = "/////"

                    if hatch is not None:
                        zorder = zorder + 2

                    bar = ax.barh(
                        i,
                        bar_width,
                        fc=model_color[model],
                        edgecolor="blue9",
                        label=models_legend[j],
                        left=starts,
                        hatch=hatch,
                        width=1.6,
                        alpha=1,
                        zorder=zorder,
                    )

    patchs = [
        mpatches.Patch(facecolor=colors_model[i], label=models_legend[i], edgecolor="blue9")
        for i in range(len(models))
    ]
    patchs.append(
        mpatches.Patch(
            facecolor="none", label="overlap", hatch="/////", edgecolor="blue9"
        )
    )
    fig.legend(patchs, ncols=3, loc="b", frame=False)
    return fig


# %%
def bar_side(MMLEA_slope, xlim=(-0.3, 0.3)):
    gs = pplt.GridSpec(nrows=1, ncols=2)
    fig = pplt.figure(refheight=4.2, refwidth=2.2, span=False, share="labels")

    models = MMLEA_slope.model.values
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
        "Dec",
        "Jan",
        "Feb",
        "Mar",
        "Apr",
        "May",
        "Jun",
        "Jul",
        "Aug",
        "Sep",
        "Oct",
        "Nov",
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
        mpatches.Patch(facecolor=colors_model[i], label=models_legend[i], edgecolor="blue9")
        for i in range(len(models))
    ]
    patchs.append(
        mpatches.Patch(
            facecolor="none", label="overlap", hatch="/////", edgecolor="blue9"
        )
    )
    fig.legend(patchs, ncols=3, loc="b", frame=False)
    return fig


#%%
def bar_dual(MMLEA_slope, xlim=(-0.35, 0.65)):

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
        "Dec",
        "Jan",
        "Feb",
        "Mar",
        "Apr",
        "May",
        "Jun",
        "Jul",
        "Aug",
        "Sep",
        "Oct",
        "Nov",
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
                yticks=np.arange(12),
                yticklabels=months,
            )

            for i, month in enumerate(months):
                data_month = MMLEA_slope.sel(
                    mode=mode, extr_type=extr_type, month=month
                )
                data_month = data_month.sortby(data_month)
                data_cum_pos = 0
                data_cum_neg = 0

                # if there are more than 4 elements in data_month whose values are in the same sign, hatch
                
                if len(data_month.where(data_month > 0).dropna("model")) == 4:
                    hatches = [None, '/////', '/////', '/////', '/////']
                elif len(data_month.where(data_month < 0).dropna("model")) == 4:
                    hatches = ['/////', '/////', '/////', '/////', None]
                elif len(data_month.where(data_month > 0).dropna("model")) == 5:
                    hatches = ['/////', '/////', '/////', '/////', '/////']
                else:
                    hatches = [None, None, None, None, None]
                
                for j, model in enumerate(models):
                    data_model = data_month.sel(model=model).values
                    bar_width = data_model
                    if data_model > 0:
                        data_cum_pos += data_model
                        data_cum_neg += 0
                        start = data_cum_pos - bar_width

                    else:
                        data_cum_pos += 0
                        data_cum_neg += data_model
                        start = data_cum_neg - bar_width

                    h = np.where(data_month.model.values == model)[0][0]

                    bar = ax.barh(
                        i,
                        bar_width,
                        fc=model_color[model],
                        edgecolor="blue9",
                        label=models_legend[j],
                        left=start,
                        width=1.6,
                        alpha=1,
                        hatch=hatches[h],
                    )

            # vline over x = 0
            ax.axvline(0, color="k", linewidth=0.4, linestyle="--")

    patchs = [
        mpatches.Patch(facecolor=colors_model[i], label=models_legend[i], edgecolor="blue9")
        for i in range(len(models))
    ]
    patchs.append(
        mpatches.Patch(
            facecolor="none", label="4+ models agree", hatch="/////", edgecolor="blue9"
        )
    )
    fig.legend(patchs, ncols=3, loc="b", frame=False)
    return fig


# %%
