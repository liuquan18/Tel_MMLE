# %%
import xarray as xr
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
import proplot as pplt
from matplotlib.lines import Line2D

from matplotlib.ticker import MaxNLocator
import matplotlib.patches as mpatches


import src.plots.composite_plot as composite_plot
import src.plots.extreme_plot as extplt
import src.plots.statistical_overview as stat_overview
import src.obs.era5_extreme_change as era5_extreme_change
import src.extreme.extreme_count_troposphere as ext_profile
#%%

CMIP6_first_profile, CMIP6_last_profile = ext_profile.extreme_count_MPI('MPI_GE_CMIP6')
CR20_plevs = [1000, 975, 950, 925, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200]
CR20_mean_first_profile, CR20_mean_last_profile = ext_profile.extreme_count_rean('CR20', CR20_plevs)

# %%
fig, axes = pplt.subplots(
    space=0,
    width=150 / 25.4,
    height=150 / 25.4,
    hspace=5,
    wspace=2,
    nrows=2,
    ncols=2,
    sharex=False,
    sharey=False,
)
axes.format(
    abc=True,
    abcloc="ul",
    abcstyle="a",
    xlabel="Extreme event count",
    ylabel="Pressure level (hPa)",
    ylim=(1000, 200),
    xminorticks="null",
    yminorticks="null",
    grid=False,
    toplabels=("pos", "neg"),
    # leftlabels=("NAO", "EA"),
    # xlocator=20,
)

ens_labels = ["first10", "last10"]
for i, count in enumerate([CMIP6_first_profile, CMIP6_last_profile]):
    extplt.plot_extreme_count(
        count.sel(mode = 'NAO', extr_type ="pos"),
        axes[0, 0],
        label=ens_labels[i],
    )
    extplt.plot_extreme_count(
        count.sel(mode = 'NAO', extr_type ="neg"),
        axes[0, 1],
        label=ens_labels[i],
    )


mean_labels = ["first40", "last40"]
for i, count in enumerate([CR20_mean_first_profile, CR20_mean_last_profile]):
    extplt.plot_extreme_count(
        count.sel(mode = 'NAO', extr_type ="pos"),
        axes[1, 0],
        label=mean_labels[i],
    )
    extplt.plot_extreme_count(
        count.sel(mode = 'NAO', extr_type ="neg"),
        axes[1, 1],
        label=mean_labels[i],
    )
axes[0,0].set_xlim(0,2.5)
axes[0,1].set_xlim(0,2.5)

axes[1,0].set_xlim(0,7)
axes[1,1].set_xlim(0,7)

# axes[0,0].set_xticks(np.arange(1,4.1,0.5))
# axes[0,1].set_xticks(np.arange(1,4.1,0.5))

# no yticks ylabel for the second columns
for ax in axes[:, 1]:
    ax.format(ylabel="")
plt.setp(axes[0,1].get_yticklabels(), visible=False)
plt.setp(axes[1,1].get_yticklabels(), visible=False)

axes[0, 0].legend(loc="lr", ncols=1, frame=True)
axes[1, 0].legend(loc="lr", ncols=1, frame=True)
for ax in axes:
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)


# %%
