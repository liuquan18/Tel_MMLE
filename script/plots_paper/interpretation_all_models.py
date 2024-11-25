# %%
import xarray as xr
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import os

from matplotlib.lines import Line2D
from matplotlib.patches import Patch
# %%
def read_corresponding(model):
    NAO_pos_jet_GB_correspond = pd.read_csv(f'/work/mh0033/m300883/Tel_MMLE/data/{model}/mechnisms/NAO_pos_jet_GB_correspond.csv', index_col = 0)
    NAO_neg_jet_GB_correspond = pd.read_csv(f'/work/mh0033/m300883/Tel_MMLE/data/{model}/mechnisms/NAO_neg_jet_GB_correspond.csv', index_col = 0)

    
    # decade equals mimimum
    NAO_pos_jet_GB_correspond_first = NAO_pos_jet_GB_correspond[NAO_pos_jet_GB_correspond['decade'] == NAO_pos_jet_GB_correspond['decade'].min()]
    NAO_neg_jet_GB_correspond_first = NAO_neg_jet_GB_correspond[NAO_neg_jet_GB_correspond['decade'] == NAO_neg_jet_GB_correspond['decade'].min()]

    NAO_pos_jet_GB_correspond_last = NAO_pos_jet_GB_correspond[NAO_pos_jet_GB_correspond['decade'] == NAO_pos_jet_GB_correspond['decade'].max()]
    NAO_neg_jet_GB_correspond_last = NAO_neg_jet_GB_correspond[NAO_neg_jet_GB_correspond['decade'] == NAO_neg_jet_GB_correspond['decade'].max()]

    
    NAO_pos_jet_GB_correspond_first['period'] = 'first'
    NAO_neg_jet_GB_correspond_first['period'] = 'first'

    NAO_pos_jet_GB_correspond_last['period'] = 'last'
    NAO_neg_jet_GB_correspond_last['period'] = 'last'

    
    NAO_pos_corresponding = pd.concat([NAO_pos_jet_GB_correspond_first, NAO_pos_jet_GB_correspond_last])
    NAO_neg_corresponding = pd.concat([NAO_neg_jet_GB_correspond_first, NAO_neg_jet_GB_correspond_last])

    return NAO_pos_corresponding, NAO_neg_corresponding


# %%
def read_dec(model):
    
    NAO_pos_df_dec = pd.read_csv(f'/work/mh0033/m300883/Tel_MMLE/data/{model}/mechnisms/NAO_pos_df.csv', index_col = 0)
    NAO_neg_df_dec = pd.read_csv(f'/work/mh0033/m300883/Tel_MMLE/data/{model}/mechnisms/NAO_neg_df.csv', index_col = 0)

    
    NAO_dec = pd.concat([NAO_pos_df_dec, NAO_neg_df_dec])
    return NAO_dec

# %%
def read_clim(model):
    dir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/mechnisms/"

    jet_loc_clim = xr.open_dataset(dir + 'jet_loc_clim.nc').jet_loc
    GB_clim = xr.open_dataset(dir + 'GB_clim.nc').GB
    return jet_loc_clim, GB_clim
# %%
corresponding_dfs = {}
dec_dfs = {}
jet_loc_clims = {}
GB_clims = {}
# %%
ens_size = [100, 50, 40, 30, 20]
for i, model in enumerate(["MPI_GE", "CanESM2", "CESM1_CAM5", "GFDL_CM3", "MK36"]):
    NAO_pos_corresponding, NAO_neg_corresponding = read_corresponding(model)
    NAO_dec = read_dec(model)
    jet_loc_clim, GB_clim = read_clim(model)
    
    corresponding_dfs[model] = (NAO_pos_corresponding, NAO_neg_corresponding)
    jet_loc_clims[model] = jet_loc_clim
    GB_clims[model] = GB_clim

    dec_dfs[model] = NAO_dec
    # divide by ensemble size
    NAO_dec['jet_loc_north'] = NAO_dec['jet_loc_north'] / ens_size[i]
    NAO_dec['GB_above'] = NAO_dec['GB_above'] / ens_size[i]
    NAO_dec['extreme_count'] = NAO_dec['extreme_count'] / ens_size[i]

# %%
fig, axes = plt.subplots(5, 2, figsize=(8, 10))

for i, model in enumerate(["MPI_GE", "CanESM2", "CESM1_CAM5", "GFDL_CM3", "MK36"]):
    NAO_pos_corresponding, NAO_neg_corresponding = corresponding_dfs[model]
    NAO_dec = dec_dfs[model]
    jet_loc_clim = jet_loc_clims[model]
    GB_clim = GB_clims[model]
    
    ax = axes[i, 0]
    sns.kdeplot(
        data=NAO_pos_corresponding,
        x="jet_loc",
        y="GB",
        ax=ax,
        hue="period",
        fill=True,
        alpha=0.5,
        common_norm=True,
        legend=False,
    )


    sns.kdeplot(
        data=NAO_neg_corresponding,
        x="jet_loc",
        y="GB",
        ax=ax,
        hue="period",
        fill=False,
        alpha=0.5,
        common_norm=True,
        legend=False,
    )

    # vline and hline


    # vline for jet_clim_first10 and jet_clim_last10
    ax.axvline(
        jet_loc_clim.isel(time = 0), linestyle="--", label="jet climatology (1850-1859)", color="C0"
    )
    ax.axvline(
        jet_loc_clim.isel(time = -1), linestyle="--", label="jet climatology (2090-2099)", color="C1"
    )
    # hline for GB_clim_first10 and GB_clim_last10
    ax.axhline(
        GB_clim.isel(time = 0), linestyle="--", label="GB climatology (1850-1859)", color="C0"
    )
    ax.axhline(
        GB_clim.isel(time = -1), linestyle="--", label="GB climatology (2090-2099)", color="C1"
    )


    ax.set_xlim(35, 65)

    ax.set_title('')
    ax.set_xlabel("jet_loc")
    ax.set_ylabel(model)


    ax2 = axes[i, 1]
    sns.scatterplot(
        data = NAO_dec,
        x = 'jet_loc_north',
        y = 'GB_above',
        hue = 'decade',
        style = 'phase',
        legend = False,
        size = 'extreme_count',
        # sizes = (10, 500),
        ax = ax2
    )

    ax2.set_xlim(-0.4, 3.4)
    ax2.set_ylim(-0.4, 3.4)
    # a b c d
    ax2.text(-0.7, 3.2, f"{chr(102 + i)}", fontsize=12, fontweight='bold')
    ax.text(32, 5.74, f"{chr(97 + i)}", fontsize=12, fontweight='bold')
    ax2.set_ylabel('')

axes[-1, 0].set_xlabel("eddy-driven jet location ($\degree$N)")
axes[-1, 1].set_xlabel("count of jet stream location (> 1.5 std)")
axes[2,1].set_ylabel("count of GB index (> 1.5 std)")

plt.savefig('/work/mh0033/m300883/Tel_MMLE/docs/source/plots/mechism/all_models_jet_GB.pdf')
    # plt.tight_layout(h_pad= -0.6)
# %%
