#%%
import matplotlib.pyplot as plt
import proplot as pplt
import numpy as np
import scipy.ndimage as ndimage


# %%
def NAO_EA_hist2d(pcs, bins =20,levels = np.arange(0,0.31,0.03),smooth = 0):
    """pc_all, pc_first,pc_last"""

    pc_all = pcs.stack(com=("time", "ens")).squeeze()
    pc_first = pcs.isel(time=slice(0, 10)).stack(com=("time", "ens")).squeeze()
    pc_last = pcs.isel(time=slice(-10, None)).stack(com=("time", "ens")).squeeze()

    fig = pplt.figure(space=0, refwidth="25em", wspace=3, hspace=3)
    fig.format(
        abc=True,
        abcloc="ul",
        abcstyle="a",
        suptitle="hist2d of NAO and EA",
    )

    axes = fig.subplots(nrows=1, ncols=3)
    for i, pc in enumerate([pc_all,pc_first, pc_last]):
        h,x,y = np.histogram2d(pc.sel(mode="NAO"), pc.sel(mode="EA"), bins=20,density=True)

        # smooth
        h = ndimage.gaussian_filter(h, sigma=smooth)

        m = axes[i].contourf(x[:-1], y[:-1], h.T, levels = levels,extend = 'both')
        axes[i].format(
            xlabel="NAO",
            ylabel="EA",
            grid=False,
            title = 'all' if i == 0 else 'first' if i == 1 else 'last'
        )
        axes[i].format(
            xlim = (-3.5,3.5),
            ylim = (-3.5,3.5),
        )
        # add colorbar
    fig.colorbar(m,         
        loc="b",
        label="count",
        length=0.8,
        width=0.1,
        shrink=0.8,)
    