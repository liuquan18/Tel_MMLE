import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import proplot as pplt
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches


def add_legend(fig):
    """
    add fake legend.
    """
    patch = mpatches.Patch(color='grey', label='first10 years')   
    line = Line2D([0], [0], label='last10 years', color='k')
    handles = [patch,line]
    fig.legend(handles,ncols = 1,title = 'periods')

def plot_hist(first, last, std_type,hlayers,bins = np.arange(-4,4.1,1)):
    fig = pplt.figure(space = 0, refwidth = "20em")
    axes = fig.subplots(nrows = 1, ncols = 2)
    axes.format(
        abc = 'a',
        abcloc = 'ul',
        xticks = 5,
        xlocator = np.arange(-4,4.1,2),
        xminorlocator = 'null',
        yminorlocator = 'null',
        suptitle = f"PDF of index at {hlayers/100:.0f}hpa standardizised with {std_type}"
    )

    modes = ['NAO','EA']
    for i, ax in enumerate(axes):

        first_hist = first.sel(mode = modes[i],hlayers = hlayers).plot.hist(bins = bins,
                                histtype = 'stepfilled',
                                color = 'grey',
                                legend_kw = {'order':'F','title':'first10'},
                                ax = ax)
        last_hist = last.sel(mode = modes[i],hlayers = hlayers).plot.hist(bins =bins,
                                histtype = 'step',
                                color = 'k',
                                linewidth = 1.,
                                legend_kw = {'order':'F','title':'last10'},
                                ax = ax)
        ims = [first_hist, last_hist]

    axes[0].format(title = 'NAO',xlabel = 'std')
    axes[1].format(title = 'EA',xlabel = 'std')
    for ax in axes:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    add_legend(fig)