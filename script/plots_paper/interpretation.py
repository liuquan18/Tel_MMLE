# %%
import xarray as xr
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import os
from src.mechnisms.mechisms import process_GB, process_jet, correlate_jet_GB

from matplotlib.lines import Line2D
from matplotlib.patches import Patch
# %%
model = 'MPI_GE'
# %%
jet_loc, jet_loc_clim, jet_loc_std, jet_loc_north, jet_north_decade = process_jet(model)
GB, GB_clim, GB_std,GB_above, GB_above_decade = process_GB(model)

#%%
NAO_pos_jet_GB_correspond = pd.read_csv(f'/work/mh0033/m300883/Tel_MMLE/data/{model}/mechnisms/NAO_pos_jet_GB_correspond.csv', index_col = 0)
NAO_neg_jet_GB_correspond = pd.read_csv(f'/work/mh0033/m300883/Tel_MMLE/data/{model}/mechnisms/NAO_neg_jet_GB_correspond.csv', index_col = 0)

#%%
# decade equals mimimum
NAO_pos_jet_GB_correspond_first = NAO_pos_jet_GB_correspond[NAO_pos_jet_GB_correspond['decade'] == NAO_pos_jet_GB_correspond['decade'].min()]
NAO_neg_jet_GB_correspond_first = NAO_neg_jet_GB_correspond[NAO_neg_jet_GB_correspond['decade'] == NAO_neg_jet_GB_correspond['decade'].min()]

NAO_pos_jet_GB_correspond_last = NAO_pos_jet_GB_correspond[NAO_pos_jet_GB_correspond['decade'] == NAO_pos_jet_GB_correspond['decade'].max()]
NAO_neg_jet_GB_correspond_last = NAO_neg_jet_GB_correspond[NAO_neg_jet_GB_correspond['decade'] == NAO_neg_jet_GB_correspond['decade'].max()]

#%%
NAO_pos_jet_GB_correspond_first['period'] = 'first'
NAO_neg_jet_GB_correspond_first['period'] = 'first'

NAO_pos_jet_GB_correspond_last['period'] = 'last'
NAO_neg_jet_GB_correspond_last['period'] = 'last'

#%%
NAO_pos_corresponding = pd.concat([NAO_pos_jet_GB_correspond_first, NAO_pos_jet_GB_correspond_last])
NAO_neg_corresponding = pd.concat([NAO_neg_jet_GB_correspond_first, NAO_neg_jet_GB_correspond_last])


#%%
jet_loc_clim_first = jet_loc_clim.isel(time = 0)
GB_clim_first = GB_clim.isel(time = 0)

jet_loc_clim_last = jet_loc_clim.isel(time = -1)
GB_clim_last = GB_clim.isel(time = -1)

# %%
NAO_pos_df_dec = pd.read_csv(f'/work/mh0033/m300883/Tel_MMLE/data/{model}/mechnisms/NAO_pos_df.csv', index_col = 0)
NAO_neg_df_dec = pd.read_csv(f'/work/mh0033/m300883/Tel_MMLE/data/{model}/mechnisms/NAO_neg_df.csv', index_col = 0)

#%%
NAO_dec = pd.concat([NAO_pos_df_dec, NAO_neg_df_dec])
# %%
jet_GB_corr = correlate_jet_GB(jet_loc, GB)
# %%
fig = plt.figure(figsize=(10, 10))
gs = fig.add_gridspec(6,2)

ax1 = fig.add_subplot(gs[0:2, 0])
ax2 = fig.add_subplot(gs[2:4, 0])
ax3 = fig.add_subplot(gs[4:6, 0])
ax4 = fig.add_subplot(gs[0:3, 1])
ax5 = fig.add_subplot(gs[3:6, 1])


jet_loc_clim.plot(ax=ax1, label='Climatology', x='time', add_legend=False, color='black')
# clean the tile of ax
ax1.set_title('')

ax1.fill_between(jet_loc_clim.time, jet_loc_clim - jet_loc_std, jet_loc_clim + jet_loc_std, color='gray', alpha=0.5, label='std')
ax1.set_ylim(42, 57)
ax1.grid(False)  # Remove gridlines
ax_twin1 = ax1.twinx()
jet_north_decade.drop_vars('lon').plot(ax=ax_twin1, color='red', label='count of jet Norther than 1.5 std')
ax_twin1.grid(False)  # Remove gridlines

lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax_twin1.get_legend_handles_labels()
ax1.legend(lines + lines2, labels + labels2, loc='upper left', frameon=False)

ax1.set_ylabel('Latitude')
ax_twin1.set_ylabel('count')
ax1.set_title('Eddy-driven jet stream location')

GB_clim.plot(ax=ax2, label='Climatology', x='time', add_legend=False, color='black')
# clean the tile of ax
ax2.set_title('')

ax2.fill_between(GB_clim.time, GB_clim - GB_std, GB_clim + GB_std, color='gray', alpha=0.5, label='std')
ax2.set_ylim(5.45, 5.7)
ax2.grid(False)  # Remove gridlines
ax_twin2 = ax2.twinx()
GB_above_decade.plot(ax=ax_twin2, color='red', label='count of GB above 1.5 std')
ax_twin2.grid(False)  # Remove gridlines

lines, labels = ax2.get_legend_handles_labels()
lines2, labels2 = ax_twin2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc='upper left', frameon=False)

ax2.set_ylabel('GB proxy (km)')
ax_twin2.set_ylabel('count')
ax2.set_title('Greenland Blocking Index')

jet_GB_corr.plot(ax=ax3, color = 'k')
ax3.set_ylabel("Correlation")
ax3.set_title("Jet Stream location and Blocking correlation")
# without grids
ax3.grid(False)


####### ax4
sns.kdeplot(
    data=NAO_pos_corresponding,
    x="jet_loc",
    y="GB",
    ax=ax4,
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
    ax=ax4,
    hue="period",
    fill=False,
    alpha=0.7,
    common_norm=True,
    legend=False,
)



# vline for jet_clim_first10 and jet_clim_last10
ax4.axvline(
    jet_loc_clim_first, linestyle="--", label="jet climatology (1850-1859)", color="C0"
)
ax4.axvline(
    jet_loc_clim_last, linestyle="--", label="jet climatology (2090-2099)", color="C1"
)
# hline for GB_clim_first10 and GB_clim_last10
ax4.axhline(
    GB_clim_first, linestyle="--", label="GB climatology (1850-1859)", color="C0"
)
ax4.axhline(
    GB_clim_last, linestyle="--", label="GB climatology (2090-2099)", color="C1"
)

ax4.set_xlabel(r"Eddy-driven jet stream location ($\degree$N)")
ax4.set_ylabel("Greenland Blocking Index proxy (km)")

# ax4.set_xlim(35, 65)
ax4.set_ylim(5.35, 5.8)

# Create custom legend
legend_elements = [
    Line2D([0], [0], color="C0", lw=2, label="first10 (1850-1859)", linestyle="-"),
    Line2D([0], [0], color="C1", lw=2, label="last10 (2090-2099)", linestyle="-"),
    Patch(facecolor="C0", edgecolor="k", label="NAO (positive)", alpha=0.5),
    Patch(facecolor="none", edgecolor="k", label="NAO (negative)", alpha=0.5),
]


ax4.legend(handles=legend_elements, loc="upper right", frameon=False)

####### ax5 

sns.scatterplot(
    data = NAO_dec,
    x = 'jet_loc_north',
    y = 'GB_above',
    hue = 'decade',
    style = 'phase',
    legend = 'brief',
    size = 'extreme_count',
    sizes = (10, 500),
    ax = ax5
)
ax5.set_xlabel("count of jet stream location (> 1.5 std)")
ax5.set_ylabel("count of Greenland blocking (> 1.5 std)")


ax5.legend(loc='upper right', frameon=False, ncol = 2)


plt.tight_layout()

# label a, b, c, d, e
ax1.text(-0.1, 1.1, 'a', transform=ax1.transAxes,
        fontsize=14, fontweight='bold', va='top', ha='right')
ax2.text(-0.1, 1.1, 'b', transform=ax2.transAxes,
        fontsize=14, fontweight='bold', va='top', ha='right')
ax3.text(-0.1, 1.1, 'c', transform=ax3.transAxes,
        fontsize=14, fontweight='bold', va='top', ha='right')
ax4.text(-0.1, 1.1, 'd', transform=ax4.transAxes,
        fontsize=14, fontweight='bold', va='top', ha='right')
ax5.text(-0.1, 1.1, 'e', transform=ax5.transAxes,
        fontsize=14, fontweight='bold', va='top', ha='right')



plt.savefig(f"//work/mh0033/m300883/Tel_MMLE/docs/source/plots/mechism/{model}_NAO_jet_GB.pdf")

# %%
