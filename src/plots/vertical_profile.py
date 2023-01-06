import xarray as xr
import matplotlib.pyplot as plt
import proplot as pplt


def plot_period(exts, ax):
    """plot vertical profile of two peiods"""
    y = (exts.hlayers.values) / 100
    styles = ['-','--','dotted']
    labels = exts.compare.values

    for i,x_period in enumerate(exts):
        x = x_period.values
        ax.plot(x,y,linestyle = styles[i],color = 'k', label = labels[i])



def plot_diff(diff_ec, ax):
    """
    plot the difference between first10 and last10 years.
    """
    x_diff_pos = diff_ec.sel(extr_type="pos").values
    x_diff_neg = diff_ec.sel(extr_type="neg").values
    y = (diff_ec.hlayers.values) / 100
    diff_pos_profile = ax.plot(
        x_diff_pos, y, linestyle="-", color="k", label="Positive"
    )
    diff_neg_profile = ax.plot(
        x_diff_neg, y, linestyle="--", color="k", label="Negative"
    )


def plot_vertical_profile(ext_count, mode,std_type):
    """
    plot the vertical profile of extreme counts
    **Arguments**
        *ecs* the extreme counts
    """
    ext_count = ext_count.sel(mode = mode)
    ext_count_diff = ext_count.isel(compare = -1) - ext_count.isel(compare = 0) 

    fig = pplt.figure(space=0, refwidth="20em")
    axes = fig.subplots(nrows=1, ncols=3)
    axes.format(
        abc="a",
        abcloc="ul",
        xlocator=[-20, -10, 0, 10, 20, 30, 40],
        xminorticks="null",
        yminorticks="null",
        suptitle=f"{mode} extreme event counts ",
        grid = False,
    )

    titles = ["pos", "neg", "diff"]

    # first 2 cols
    extr_types = ["pos", "neg"]
    for i, ax in enumerate(axes[:2]):
        ext_type = ext_count.sel(extr_type = extr_types[i])
        plot_period(ext_type, ax=ax)
        ax.format(title=titles[i], xlabel="extremes count", ylabel="gph/hpa")

    # the last col
    plot_diff(ext_count_diff, ax=axes[2])
    axes[2].format(title=titles[2])

    for ax in axes[:2]:
        ax.set_ylim(1000, 200)
        ax.set_xlim(5, 45)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
    axes[0].legend(loc="lr", ncols=1, title="period")
    axes[1].legend(loc="ll", ncols=1, title="period")
    axes[2].legend(loc="cr", ncols=1, title="extreme type")
    axes[2].set_ylim(1000, 200)
    axes[2].set_xlim(-10, 30)
    axes[2].spines["right"].set_visible(False)
    axes[2].spines["top"].set_visible(False)


