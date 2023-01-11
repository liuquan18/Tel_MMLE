import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns
import xarray as xr
from cartopy.mpl.geoaxes import GeoAxes
from cartopy.mpl.ticker import LatitudeFormatter, LatitudeLocator, LongitudeFormatter
from cartopy.util import add_cyclic_point
from matplotlib.colorbar import Colorbar
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import AxesGrid, make_axes_locatable

from cartopy.mpl.geoaxes import GeoAxes
from cartopy.mpl.ticker import LatitudeFormatter, LatitudeLocator, LongitudeFormatter
from cartopy.util import add_cyclic_point
from matplotlib.colorbar import Colorbar
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import AxesGrid, make_axes_locatable
from cartopy.mpl.geoaxes import GeoAxes
from cartopy.mpl.ticker import LatitudeFormatter, LatitudeLocator, LongitudeFormatter
from cartopy.util import add_cyclic_point
from matplotlib.colorbar import Colorbar
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import AxesGrid, make_axes_locatable


import iris
import iris.plot as iplt
import iris.quickplot as qplt


# function to erase the white line
def erase_white_line(data):
    """
    erase the white line aroung 180 degree.
    """
    data = data.transpose(..., "lon")  # make the lon as the last dim
    dims = data.dims  # all the dims
    res_dims = tuple(dim for dim in dims if dim != "lon")  # dims apart from lon
    res_coords = [data.coords[dim] for dim in res_dims]  # get the coords

    # add one more longitude to the data
    data_value, lons = add_cyclic_point(data, coord=data.lon, axis=-1)

    # make the lons as index
    lon_dim = xr.IndexVariable(
        "lon", lons, attrs={"standard_name": "longitude", "units": "degrees_east"}
    )

    # the new coords with changed lon
    new_coords = res_coords + [lon_dim]  # changed lon but new coords

    new_data = xr.DataArray(data_value, coords=new_coords, name=data.name)

    return new_data


def buildax(ax, zorder=50, alpha_coast=0.7,alpha_grid = 0.5):
    """
    add grid coastline and gridlines
    """
    ax.set_global()
    ax.coastlines(linewidth=0.5, alpha=alpha_coast)
    gl = ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=False,
        linewidth=0.5,
        zorder=zorder,
        alpha=alpha_grid,
    )
    gl.xformatter = LongitudeFormatter(zero_direction_label=False)
    gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, 45))

    gl.ylocator = mticker.FixedLocator([20, 40, 60])
    gl.yformatter = LatitudeFormatter()


def remove_cb(ax):
    """
    remove the colorbar of the object
    """
    cb = ax.colorbar
    cb.remove()


def add_cb(
    fig, im, loc=[0.85, 0.2, 0.03, 0.6], label="label", orientation="horizontal"
):
    """
    add colobar to a fig based on the image at loc.
    """

    cbar_ax = fig.add_axes(loc)
    fig.colorbar(im, cax=cbar_ax, label=label, orientation=orientation)
