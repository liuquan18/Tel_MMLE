#%%
import matplotlib.pyplot as plt
import proplot as pplt
# %%
def NAO_EA_hist2d(pcs,bins = 30):
    """pc_all, pc_first,pc_last"""
    fig = pplt.figure(space=0, refwidth="25em", wspace=3, hspace=3)
    fig.format(
        abc=True,
        abcloc="ul",
        abcstyle="a",
        title="hist2d of NAO and EA",
    )

    axes = fig.subplots(nrows=3, ncols=1,share=False)

    pc_all = pcs.stack(com = ('time','ens')).squeeze()
    pc_fist = pcs.isel(time=(0,10)).stack(com = ('time','ens')).squeeze()
    pc_last = pcs.isel(time=(-10,None)).stack(com = ('time','ens')).squeeze()

    for i, pc in enumerate([pc_all,pc_fist,pc_last]):
        axes[i].hist2d(pc.sel(mode='NAO'),pc.sel(mode='EA'),bins=bins)
        axes[i].format(
            xlabel='NAO',
            ylabel='EA',
            grid=False,
            suptitle='hist2d of NAO and EA',
        )
