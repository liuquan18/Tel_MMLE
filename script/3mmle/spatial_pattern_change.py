#%%
import xarray as xr
import numpy as np
import pandas as pd

#%%
import src.Teleconnection.spatial_pattern as ssp
import src.Teleconnection.tools as tools
import proplot as pplt
import cartopy.crs as ccrs


# %%
data = xr.open_mfdataset(
    "/work/mh0033/m300883/Tel_MMLE/data/test_CESM/zg_processed/*.nc",
    combine="nested",
    concat_dim="ens",
    join="override",
)

#%%
data = data.rename({'plev':'hlayers'}).zg
data = data.sel(hlayers = slice(100000,20000))

periods = [slice("1920","1930"),slice("2000","2010"),slice("2090","2100")]

# %%
def spatial_pattern_change(data, periods,dim = 'hlayers'):
    """
    get the spatial pattern of the data in periods
    """
    if data.hlayers.size>1:
        data = tools.standardize(data)
    EOFs = []
    period_index = xr.IndexVariable(dims="period", data=periods)
    for period in periods:
        data_p = data.sel(time=period)
        EOF = vertical_spatial_pattern(data_p,dim = dim)
        EOFs.append(EOF)

    EOFs = xr.concat(EOFs, dim=period_index)
    EOFs['period'] = ['0C','2C','4C']

    return EOFs


def vertical_spatial_pattern(data,dim = 'hlayers'):
    """
    get the spatial pattern for all levels.
    """
    data = data.stack(com = ("time","ens"))
    eofs = data.groupby(dim,squeeze = True).apply(spatial_pattern)
    return eofs


def spatial_pattern(data):
    """
    get the spatial pattern of single data
    """
    eof, _, _ = ssp.doeof(data, nmode=2, dim="com", standard=True)

    return eof

# %%

def spatial_stat(eof,mode,dim = 'lon'):
    """
    zonally mean ('lat' as the final dim) or meridional mean ('lon' as the final dim)
    **Arguments**
        *eof* the eof
        *mode* the mode 
        *dim* the dim that left
    """
    if dim == 'lon':
        bins = np.arange(-90, 41, 5) # lon bins
        labels = np.arange(-90,36,5)  # lon label
        average_dim = 'lat'
    elif dim == 'lat':
        bins = np.arange(20,81,4)
        labels = np.arange(22,81,4)
        average_dim = 'lon'

    data = eof.sel(mode = mode)
    data_height = (data.groupby_bins(dim,bins = bins, labels = labels)
                    .mean(dim = average_dim )
    )
    return data_height
    
# %%


# PLOT maps
def spatial_pattern_maps(eofs, hlayers = 50000,levels=np.arange(-2.0, 2.1, 0.4)):
    """
    rows as modes
    cols in different periods
    """
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
        toplabels=("0C","2C","4C"),
        leftlabels=("NAO", "EA"),
        suptitle=f"spatial patterns on {hlayers/100:.0f}hpa",
    )

    for r, mode in enumerate(eofs.mode):
        for c, period in enumerate(eofs.period):

            eof = eofs.sel(hlayers = hlayers, mode = mode, period = period)
            map = axes[r,c].contourf(
                eof,
                x = 'lon',
                y = 'lat',
                levels = levels,
                extend = 'both',
                transform = ccrs.PlateCarree(),
                cmap = "RdBu_r"
            )

    fig.colorbar(map, loc = "r", pad = 3, title = "gph/std")
# %%
