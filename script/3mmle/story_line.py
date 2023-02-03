# %%
# import xarray, numpy, pandas, proplot
import xarray as xr
import numpy as np
import pandas as pd
import proplot as pplt
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import src.plots.plot_violin as violin_plots
import seaborn as sns


# extremes
import src.extreme.extreme_ci as extreme

# warming stage
import src.warming_stage.warming_stage as warming_stage

#%%
# load data
odir = "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct/"
pc_dir = odir + "EOF_result/ind_first_pc.nc"
eof_dir = odir + "EOF_result/ind_first_eof.nc"
fra_dir = odir + "EOF_result/ind_first_fra.nc"

tsurf_dir = odir + "ts_processed/tsurf_mean.nc"

pc = xr.open_dataset(pc_dir).pc
eof = xr.open_dataset(eof_dir).eof
fra = xr.open_dataset(fra_dir).exp_var
tsurf = xr.open_dataset(tsurf_dir).tsurf
# %%
# split into first 10 and last 10 years
periods_pc,periods = warming_stage.split_period(pc, compare = 'CO2')
first_pc, last_pc = periods_pc[0], periods_pc[1]

#%%
# extreme events count
first_count = extreme.extreme_count_xr(first_pc)
last_count = extreme.extreme_count_xr(last_pc)


#*************** Figs ************************
#%%
# Fig 1 spatial patterns and statistics of the pcs
# rows for 'NAO' and 'EA'
# colums for 'spatial map','pc hist','violion vertical profile'.

fig = pplt.figure(space=0, refwidth="25em", wspace=3, hspace=3)
fig.format(
    abc=True,
    abcloc="ul",
    abcstyle="a",
    title="spatial patterns and statistics of the pcs",
    leftlabels=("NAO", "EA"),
)


gs = pplt.GridSpec(
    ncols=3,
    nrows=2,
    wspace=2,
    wratios=(
        1,
        0.8,
        1
    ),
)
modes = ["NAO", "EA"]

for i, mode in enumerate(modes):
    # plot spatial map
    spatial_ax = fig.add_subplot(gs[i, 0], proj="ortho", proj_kw={"lon_0": -20, "lat_0": 60})
    map = spatial_ax.contourf(
        eof.sel(mode=mode, hlayers=50000), levels=np.arange(-1, 1.1, 0.1), extend="both"
    )
    spatial_ax.format(
        lonlines=20,
        latlines=30,
        coast=True,
        coastlinewidth=0.5,
        coastcolor="charcoal",
        title = mode + f"({fra.sel(mode = mode,hlayers = 50000):.0%})"
    )

    # plot pc hist
    hist_ax = fig.add_subplot(gs[i, 1])

    bins = np.arange(-4,4.1,1)

    first_hist = first_pc.sel(mode = mode,hlayers = 50000).plot.hist(bins = bins,
                            histtype = 'stepfilled',
                            color = 'grey',
                            legend_kw = {'order':'F','title':'first10'},
                            ax = hist_ax)
    last_hist = last_pc.sel(mode = mode,hlayers = 50000).plot.hist(bins =bins,
                            histtype = 'step',
                            color = 'k',
                            linewidth = 1.,
                            legend_kw = {'order':'F','title':'last10'},
                            ax = hist_ax)
    # legend
    patch = mpatches.Patch(color='grey', label='first10 years')   
    line = Line2D([0], [0], label='last10 years', color='k')
    handles = [patch,line]
    
    hist_ax.format(
        grid=False,
        yminorticks="null",
        xminorticks="null",
        title = mode
    )
    hist_ax.spines['right'].set_visible(False)
    hist_ax.spines['top'].set_visible(False)

    # plot violin
    violin_ax = fig.add_subplot(gs[i, 2])
    df = violin_plots.xr2df(first_pc, last_pc, mode=modes[i], compare='CO2')
    g = sns.violinplot(
        data=df,
        y="hlayers",
        x=modes[i],
        hue="period",
        kind="violin",
        palette="pastel",
        orient="h",
        ax=violin_ax,
        split=False,
        dodge=True,
        linewidth=1,
    )
    g.axes.legend().remove()
    violin_ax.format(
        grid=False,
        yminorticks="null",
        xminorticks="null",
        title=mode,
        xlabel = 'std',
        ylabel = 'gph/hPa',
        xlim = (-5,5)
    )
    violin_ax.spines['left'].set_visible(False)
    violin_ax.spines['right'].set_visible(False)
    violin_ax.spines['top'].set_visible(False)

    if i == 1:
        spatial_ax.colorbar(map, loc="b", title="std", ticks=0.2, pad=2)
        hist_ax.legend(handles,loc = 'b', ncols = 2,title = 'periods')
        violin_ax.legend(loc="b", ncols=2, title="periods")
    plt.savefig('/work/mh0033/m300883/Tel_MMLE/docs/source/plots/story_line/Fig1.png')
# %%
# Fig 2  extreme event count profile
