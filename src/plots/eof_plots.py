import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.path as mpath
from matplotlib.colorbar import Colorbar
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D


import cartopy.crs as ccrs
import cartopy.mpl.ticker as cticker
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.geoaxes import GeoAxes
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter,
                                LatitudeLocator)

import seaborn as sns
import numpy as np
import xarray as xr
import pandas as pd

import src.Teleconnection.spatial_pattern as ssp
import src.Teleconnection.index_statistic as sis


def exp_time_height(exp_plot):
    """
    draw contourf of explained variance as function of time and height.
    """

    fig,axes = plt.subplots(1,3,figsize = (14,4),dpi = 150)
    naoexp= exp_plot.sel(hlayers = slice(200,1000),mode = 'NAO').plot(x = 'time',y = 'hlayers',
                                                        ax = axes[0],
                                                        levels = np.arange(30,46,2),
                                                        extend = 'both'
                                                                    )
    eaexp = exp_plot.sel(hlayers = slice(200,1000),mode = 'EA').plot(x = 'time',y = 'hlayers',
                                                        ax = axes[1],
                                                        # levels = np.arange(12,20,0.5),
                                                        extend = 'both')

    axes[0].set_ylim(1000,200)
    axes[1].set_ylim(1000,200)


    axes[0].set_ylabel("gph/hpa")
    axes[1].set_ylabel(None)

    naocb = naoexp.colorbar.ax
    naocb.set_ylabel('exp/%')

    eacb = eaexp.colorbar.ax
    eacb.set_ylabel('exp/%')

    axes[0].set_title("NAO explained variance")
    axes[1].set_title("EA explaeined variance")

def axbuild(ax):
    
    theta = np.linspace(0,2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    
    
    ax.coastlines()
    gl=ax.gridlines(crs = ccrs.PlateCarree(),draw_labels=False)
    gl.xformatter = LongitudeFormatter(zero_direction_label=False)
    gl.xlocator = mticker.FixedLocator(np.arange(-180,180,45))

    gl.ylocator = mticker.FixedLocator([20,40,60])
    gl.yformatter = LatitudeFormatter()

    ax.get_extent(crs = ccrs.PlateCarree())
    ax.set_extent([-180,180,20,90],crs = ccrs.PlateCarree())
    ax.set_boundary(circle, transform=ax.transAxes)
    return ax

def visu_eofspa(eofs,plev = [50000,85000],levels = np.arange(-1,1.1,0.2)):
    """
    two heights, both NAO and EA
    """
    fig,axes = plt.subplots(2,2,figsize = (8,8),dpi = 500,
                        subplot_kw={'projection':
                                    ccrs.LambertAzimuthalEqualArea(
                                        central_longitude=0.0,
                                        central_latitude=90.0)
                                    })                     
    for ax in axes.flat:
        axbuild(ax)
    
    mode = ['NAO','EA']
        
    for i,row in enumerate(axes): # plev
        for j, col in enumerate(row):  # mode
            data = eofs.sel(hlayers = plev[i], mode = mode[j]).values
            im = col.contourf(eofs.lon,eofs.lat,data,
                        levels = levels,
                        extend = 'both',
                        transform = ccrs.PlateCarree(),
                        cmap = 'RdBu_r'
                        )
            col.set_title('plev:{:3.0f} hpa mode:{}'.format(plev[i]/100,mode[j]))
    fig.subplots_adjust(hspace = 0.05,wspace = 0.05,right = 0.8)
    cbar_ax = fig.add_axes([0.85, 0.2, 0.03, 0.6])
    fig.colorbar(im, cax=cbar_ax,label = 'eofs')

    plt.show()


def visu_eofspa_all(eofs,mode = 'EA',levels = np.arange(-1,1.1,0.2)):
    """
    one mode, all heights
    """
    cols = len(eofs.hlayers)//3+1
    fig,axes = plt.subplots(3,cols,figsize = (3*3,3*cols),
                        subplot_kw={'projection':
                                    ccrs.LambertAzimuthalEqualArea(
                                        central_longitude=0.0,
                                        central_latitude=90.0)
                                    })                     
    for i,ax in enumerate(axes.flat):
        axbuild(ax)
        if i < len(eofs.hlayers):
            data = eofs.sel(mode=mode).isel(hlayers = i).values
            im = ax.contourf(eofs.lon,eofs.lat,data,
                            levels = levels,
                            extend = 'both',
                            transform = ccrs.PlateCarree(),
                            cmap = 'RdBu_r'        
            )
            ax.set_title("eof {}".format(eofs.isel(hlayers = i).hlayers.values))
    cbar_ax = fig.add_axes([0.85, 0.2, 0.03, 0.6])
    fig.colorbar(im, cax=cbar_ax,label = 'eofs')

    
def visu_eof_single(eof,levels = np.arange(-1,1.1,0.2)):
    """
    one height,both NAO and EA
    """
    EOFmaps = eof.plot.contourf('lon','lat',col = 'mode',
                                        levels = levels,
                                        extend = 'both',
                                        subplot_kws=dict(projection = ccrs.LambertAzimuthalEqualArea(central_longitude=0.0,
                                                                                                    central_latitude=90.0),
                                                        )
                                        ,transform = ccrs.PlateCarree(),add_colorbar = True )

    for i,ax in enumerate(EOFmaps.axes.reshape(-1)):
        axbuild(ax)
        
        ax.set_title(f'mode={eof.mode[i].values}')
        
        
    fig = EOFmaps.fig
    fig.set_figheight(6)
    fig.set_figwidth(13.5)

    EOFmaps.cbar.set_label("gph/m")


def visu_spatial_type(eofs,plev,mode = 'EA'):
    """
    oen mode, three patterns
    """
    all_eof,first_eof,last_eof = [eof.sel(hlayers = plev) for eof in eofs]
    eof = xr.concat([first_eof,all_eof,last_eof],dim = 'type',coords='minimal',compat='override')
    eof['type'] = ['first','all','last']

    EOFmaps = eof.sel(mode=mode).plot.contourf('lon','lat',col = 'type',
                                        levels = np.arange(-1,1.1,0.2),
                                        extend = 'both',
                                        subplot_kws=dict(projection = ccrs.LambertAzimuthalEqualArea(central_longitude=0.0,
                                                                                                    central_latitude=90.0),
                                                        )
                                        ,transform = ccrs.PlateCarree(),add_colorbar = True )

    for i,ax in enumerate(EOFmaps.axes.reshape(-1)):
        axbuild(ax)
        
        ax.set_title(f'{eof.type[i].values}')
        
        
    fig = EOFmaps.fig
    fig.set_figheight(6)
    fig.set_figwidth(13.5)

    EOFmaps.cbar.set_label("gph/m")
    plt.suptitle(f"plev = {plev}")
    plt.show()

def visu_composite_spa(composite,plev = 50000,levels = np.arange(-2,2.1,0.4)):
    """
    two heights, both NAO and EA
    """
    fig,axes = plt.subplots(2,2,figsize = (8,8),dpi = 500,
                        subplot_kw={'projection':
                                    ccrs.LambertAzimuthalEqualArea(
                                        central_longitude=0.0,
                                        central_latitude=90.0)
                                    })                     
    for ax in axes.flat:
        axbuild(ax)
    
    mode = ['NAO','EA']
    extr_type = ['pos','neg']
        
    for i,row in enumerate(axes): # extr_type
        for j, col in enumerate(row):  # mode
            data = composite.sel(hlayers = plev,extr_type = extr_type[i], mode = mode[j]).values
            im = col.contourf(composite.lon,composite.lat,data,
                        levels = levels,
                        extend = 'both',
                        transform = ccrs.PlateCarree(),
                        cmap = 'RdBu_r'
                        )
            col.set_title(f'extreme:{extr_type[i]} mode:{mode[j]}')
    fig.subplots_adjust(hspace = 0.05,wspace = 0.05,right = 0.8)
    cbar_ax = fig.add_axes([0.85, 0.2, 0.03, 0.6])
    fig.colorbar(im, cax=cbar_ax,label = 'eofs')

    plt.show()


def tenyr_hist(data,hlayer = 50000,bins = 50):
    """
    visu 2d histplots of first-first---first-all, last-last---last-all
    """

    if hlayer == 'all':
        data = data
    else:
        data = data.loc[hlayer]
    fig,ax = plt.subplots()
    hf = sns.histplot(data = data,x = 'pc_all',y = 'pc_first',
    ax = ax,color= 'b', bins = bins,label = 'first',legend = False,alpha = 0.9)

    hl = sns.histplot(data = data,x = 'pc_all',y = 'pc_last',
    ax = ax,color = 'r', bins = bins,label = 'last',legend = False,alpha = 0.9)

    line = ax.plot(np.arange(-3,4,1),np.arange(-3,4,1),linestyle = 'dotted',color = 'k')

    blue_patch = mpatches.Patch(color='blue',label="first")
    red_patch = mpatches.Patch(color='red', label='last')

    plt.legend(handles=[blue_patch,red_patch],loc = 'upper left')


def tenyr_scatter(first,last,hlayer = 'all'):
    """
    make scaterplots of first_on_first v.s first_on_all and last_on_last v.s last_on_all.
    """

    fig, axes = plt.subplots(1,2,figsize = (8,3.5),dpi = 150)
    plt.subplots_adjust(wspace = 0.3)
    modes = ['NAO','EA']
    for i,ax in enumerate(axes):
        if hlayer=='all':
            first_data,last_data = first.loc[:,modes[i],:],last.loc[:,modes[i],:]
        else:
            first_data,last_data = first.loc[hlayer,modes[i],:],last.loc[hlayer,modes[i],:]

        scaf = sns.scatterplot(data = first_data, x = 'pc_all',y = 'pc_first',
        ax = ax,label = 'first')
        scar = sns.scatterplot(data = last_data, x = 'pc_all',y = 'pc_last',
        ax = ax,label = 'last',color = 'r',alpha=0.3)
        
        line = ax.plot(np.arange(-5,5,1),np.arange(-5,5,1),linestyle = 'dotted',color = 'k')

        ax.legend(loc = 'upper left')
        ax.set_ylabel('pc/std')
        ax.set_xlabel('pc_all/std')
    axes[0].set_title("NAO")
    axes[1].set_title("EA")
    plt.suptitle("first_first on first_all and last_last on last_all")


def tenyr_scatter_extreme(first,last,hlayer = 'all'):
    """
    make scaterplots of two first_on_first_first_on_all and last_on_last_last_on_all.
    """

    fig, axes = plt.subplots(2,2,figsize = (7,7),dpi = 150)
    plt.subplots_adjust(wspace = 0.3,hspace = 0.3)
    modes = ['NAO','EA']
    for i,row in enumerate(axes.T):
        for ax in row:
            if hlayer=='all':
                first_data,last_data = first.loc[:,modes[i],:],last.loc[:,modes[i],:]
            else:
                first_data,last_data = first.loc[hlayer,modes[i],:],last.loc[hlayer,modes[i],:]

            scaf = sns.scatterplot(data = first_data, x = 'pc_all',y = 'pc_first',
            ax = ax,label = 'first')
            scar = sns.scatterplot(data = last_data, x = 'pc_all',y = 'pc_last',
            ax = ax,label = 'last',color = 'r',alpha=0.3)
            
            line = ax.plot(np.arange(-5,5,1),np.arange(-5,5,1),linestyle = 'dotted',color = 'k')

            ax.legend(loc = 'upper left')
            ax.set_ylabel('pc/std')
            ax.set_xlabel('pc_all/std')


    axes[0,0].set_xlim(2,4)
    axes[0,0].set_ylim(2,4)

    axes[0,1].set_xlim(2,4)
    axes[0,1].set_ylim(2,4)

    axes[1,0].set_xlim(-4,-2)
    axes[1,0].set_ylim(-4,-2)
    axes[1,1].set_xlim(-4,-2)
    axes[1,1].set_ylim(-4,-2) 

    axes[0,0].set_title("NAO")
    axes[0,1].set_title("EA")


def scatter_extreme(*args,mode = 'NAO', hlayer = 'all'):
    """
    plot extreme scatter of one mode at all three periods (first10, last10, dynamic)
    **Arguments**
        *dfs* the three rows of dataframes to plot.
               for each period, dataframes of projection on first and last10 should be 
               included. e.g. [first_first_all, last_first_all]
    """
    nperiods = len(args)
    fig, axes = plt.subplots(2,nperiods,figsize = (8,5),dpi = 150) # rows for pos-neg, 
                                                                   # columns for periods.
    plt.subplots_adjust(hspace = 0.4,wspace = 0.4)
    for row in axes: # first row for positive extreme, second row for negative extreme
        for i, period in enumerate(args):
            if hlayer == 'all':
                firstPattern, lastPattern = [period_ten.loc[:,mode,:,:]
                for period_ten in period]
            else:
                firstPattern, lastPattern = [period_ten.loc[hlayer,mode,:,:]
                for period_ten in period]
            
            scatterfirst = sns.scatterplot(data = firstPattern,x = 'pc_all',y = 'pc_first',
            ax = row[i],label = 'first_pattern')
            scatterlast = sns.scatterplot(data = lastPattern,x = 'pc_all',y = 'pc_last',
            ax = row[i], label = 'last_pattern',color = 'r',alpha=0.5)

            line = row[i].plot(np.arange(-5,5,1),np.arange(-5,5,1),
            linestyle = 'dotted',color = 'k')

            row[i].legend(loc = 'upper left',fontsize = 7)
            row[i].set_ylabel("index")
            row[i].set_xlabel("index_all")
    for ax in axes[0]:
        ax.set_xlim(2,4)
        ax.set_ylim(2,4)
    for ax in axes[1]:
        ax.set_xlim(-4,-2)
        ax.set_ylim(-4,-2)
    axes[0,0].set_title("first 10 period")
    axes[0,1].set_title("last 10 period")
    axes[0,2].set_title("dynamic")

def extreme_bar(extreme_counts,mode = 'NAO',hlayer = 'all',ylim = 360):
    """
    plot the barplot of extreme counts. rows for 'pos' or 'neg'. cols for 'ind' or 'dep'
    **Arguments**
        *extreme_counts* the data for 'ind' and 'dep'.
        *mode* 'NAO' or 'EA'
    **Return**
        plots
    """
    fig,axes = plt.subplots(2,2,figsize = (6,3),dpi = 150)
    plt.subplots_adjust(hspace = 0)

    colors = ['#1f77b4', '#2ca02c', '#d62728']
    extr_type = ['pos','neg']

    for i, row in enumerate(axes):
        for j, col in enumerate(row): # ['ind' or 'dep']
            if hlayer == 'all':
                data = sis.all_layer_counts(extreme_counts[j]).loc[extr_type[i],mode]
            else:
                data = extreme_counts[j].loc[extr_type[i],mode,hlayer]
            sns.barplot(data = data, x = 'period',y = 'extreme_counts',ax = col,
            hue = 'pattern', hue_order=['first','all','last'],palette = colors)
            if i ==0:
                col.set_ylim(0,ylim)
            if i ==1:
                col.set_ylim(ylim,0)

    axes[0,0].set_ylabel("positive")
    axes[1,0].set_ylabel("negative")

    axes[0,1].set_ylabel(None)
    axes[1,1].set_ylabel(None)
    axes[0,1].get_legend().remove()
    axes[1,0].get_legend().remove()
    axes[1,1].get_legend().remove()
    axes[0,0].legend(loc = 'upper left')


def vertical_profile(extreme_counts,mode = 'NAO'):
    """
    using matplotlib to plot the vertical profile.
    solve the problem of y-axis sort.
    """

    # combine and unstack
    all = sis.combine_diff(extreme_counts,mode = mode)

    # plot  
    fig, axes = plt.subplots(1,3,figsize = (8,3),dpi = 150)
    plt.subplots_adjust(wspace = 0.3)

    colors = ['#1f77b4', '#2ca02c', '#d62728','#ff7f0e']
    patterns = ['first','all','last']
    periods =['first10','last10','diff']

    for i, ax in enumerate(axes):  # periods
        period_data = all[i]

        for j, pattern in enumerate(period_data.columns.levels[0]):
            pattern_data = period_data[pattern].sort_index()
            y = (pattern_data.index.values/100).astype(int)

            ax.plot(pattern_data['pos'], y, color = colors[j])
            ax.plot(pattern_data['neg'], y, color = colors[j],dashes = [3,3])

            ax.set_ylim(1000,200)
            if i<2:
                ax.set_xlim(0,50)
            elif i ==2:
                ax.set_xlim(-10,40)

            ax.set_title(f"{mode} {periods[i]}")

            ax.set_xlabel("extreme counts")
            if i == 0:
                ax.set_ylabel("gph/hpa")

    # legend
    custom_lines = [Line2D([0],[0],color = colors[0]),
                    Line2D([0],[0],color = colors[1]),
                    Line2D([0],[0],color = colors[2]),
                    Line2D([0],[0],color = colors[3]),
                    Line2D([0],[0],color = None,alpha = 0),
                    Line2D([0],[0],color = 'k'),
                    Line2D([0],[0],dashes = [3,3],color = 'k')]
    type_legend = axes[-1].legend(custom_lines,['first','all','last','dynamic','','pos','neg'],
    loc = 'lower right',fontsize = 6)
    axes[-1].add_artist(type_legend)

def vertical_profile_diff(ind_extre,dep_extre):
    """
    vertical profle line plots, but for difference only.
    """
    ind_diff = sis.period_diff(ind_extre).unstack([0,1])
    dep_diff = sis.period_diff(dep_extre).unstack([0,1])
    
    fig,axes = plt.subplots(2,3,figsize =(7,7),dpi = 150)
    plt.subplots_adjust(wspace=0.3,hspace=0.3)

    colors = ['#1f77b4', '#2ca02c', '#d62728','#ff7f0e']
    patterns = ['first','all','last','dynamic']
    modes = ['NAO','EA']
    inds = ['independent','dependent','difference']
    data = [ind_diff,dep_diff,ind_diff-dep_diff]

    for i, row in enumerate(axes): # modes
        for j, ax in enumerate(row):  # ind
            for l,p in enumerate(patterns):
                data_mode = data[j].xs((p,modes[i]),level = ['diff','mode'],axis = 1)\
                    .sort_index()
                y = (data_mode.index.values/100).astype(int)

                ax.plot(data_mode['pos'], y, color = colors[l])
                ax.plot(data_mode['neg'], y, color = colors[l],dashes = [3,3])

                ax.set_ylim(1000,200)
                ax.set_xlim(-15,40)
                ax.set_title(f"{modes[i]} {inds[j]}")

            if j == 0:
                axes[i,j].set_ylabel("gph/hpa")
            if i == 1:
                axes[i,j].set_xlabel("extreme counts")


    # legend
    custom_lines = [Line2D([0],[0],color = colors[0]),
                    Line2D([0],[0],color = colors[1]),
                    Line2D([0],[0],color = colors[2]),
                    Line2D([0],[0],color = colors[3]),
                    Line2D([0],[0],color = None,alpha = 0),
                    Line2D([0],[0],color = 'k'),
                    Line2D([0],[0],dashes = [3,3],color = 'k')]
    type_legend = axes[0,0].legend(custom_lines,['first','all','last','dynamic','','pos','neg'],
    loc = 'lower right',fontsize = 6)
    axes[0,0].add_artist(type_legend)

    type_legend = axes[1,0].legend(custom_lines,['first','all','last','dynamic','','pos','neg'],
    loc = 'lower right',fontsize = 6)
    axes[1,0].add_artist(type_legend)

def scatter_pattern_counts(ind_pattern,dep_pattern,mode,fit_reg = True):
    """
    scatter plot, extreme events increments v.s pattern difference
    **Arguments**
        *ind_pattern* the dataframe, with columns of extreme-event-number change 
        (compare to extreme counts with first pattern) and pattern difference from 
        first pattern.
        *dep_pattern* the same, but for dependent eof analysis.
        *mode* NAO or EA
    **Return**
        axes
    """
    fig,axes = plt.subplots(1,2,figsize = (7,4),dpi = 500)

    ind = ind_pattern.xs((slice(70000,100000),mode),level = ('hlayers','mode'))

    dep = dep_pattern.xs((slice(70000,100000),mode),level = ('hlayers','mode'))

    sns.regplot(data = ind,y = 'pos',x = 'pattern_diff',ax = axes[0],
    label = 'positive',fit_reg = fit_reg)

    sns.regplot(data = ind,y = 'neg',x = 'pattern_diff',ax = axes[0],
    label = 'negative',
    robust=True,fit_reg = fit_reg)

    sns.regplot(data = dep,y = 'pos',x = 'pattern_diff',ax = axes[1],
    label = 'positive',fit_reg = fit_reg)

    sns.regplot(data = dep,y = 'neg',x = 'pattern_diff',ax = axes[1],
    label = 'negative',fit_reg = fit_reg)



    for ax in axes:
        ax.set_xlim(0,0.22)
        ax.set_ylim(-1,20)
        ax.set_xlabel("pattern difference")
        ax.set_ylabel("extreme count increment difference")
        ax.legend(loc = 'upper left')
    axes[0].set_title("independent")
    axes[1].set_title("dependent")


def vertical_profile_slides(extreme_counts,mode = 'NAO'):
    """
    using matplotlib to plot the vertical profile.
    solve the problem of y-axis sort.
    """

    # select one single pattern for each subplot.
    all = sis.combine_diff(extreme_counts,mode = mode)
    first = all[0].xs('first',level = 'pattern',axis = 1)
    last = all[1].xs('last',level = 'pattern',axis = 1)
    dynamic = all[2].xs('dynamic',level = 'diff',axis = 1)
    data_single = pd.concat([first,last,dynamic],keys = ['first','last','dynamic'],axis = 1)

    # plot  
    fig, axes = plt.subplots(1,3,figsize = (8,3),dpi = 300)
    plt.subplots_adjust(wspace = 0.3)

    patterns = ['first','last','dynamic']
    periods =['first10','last10','diff']

    for i, ax in enumerate(axes):  # periods
        period_data = data_single[patterns[i]]
        y = (period_data.index.values/100).astype(int)
        
        ax.plot(period_data['pos'], y, c = 'k',ls = 'solid',label = 'pos')
        ax.plot(period_data['neg'], y, c = 'k',ls = 'dashed',label = 'neg')

        ax.set_ylim(1000,200)
        ax.set_xlim(-5,50)

        ax.set_title(f"{mode} {periods[i]}")

        ax.set_xlabel("extreme counts")
        if i == 0:
            ax.set_ylabel("gph/hpa")
    axes[0].legend(loc = 'lower right')
