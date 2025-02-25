#%%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
# %%
def read_eof_decade(model, fixed_pattern="decade_mpi"):
    """read eofs that is decomposed by decade"""
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/EOF_result/"
    filename = f"plev_50000_{fixed_pattern}_first_JJA_eof_result.nc"
    ds = xr.open_dataset(odir + filename)
    ds = ds.sel(mode="NAO")
    return ds
#%%
def quantile_90_10(eof, quantile=0.9):
    first_pc = eof.pc.isel(time = slice(0,30))
    # 90th quantile
    first_quantile = first_pc.quantile(quantile)

    return first_quantile.values
#%%
def percentage_le_1_5(eof, threshold = 1.5):
    first_pc = eof.pc.isel(time = slice(0, 30))
    # Calculate the fraction of values <= 1.5
    fraction_le_1_5 = (first_pc <= threshold).mean()

    # Convert fraction to percentage
    percentage_le_1_5 = fraction_le_1_5 * 100
    return percentage_le_1_5.values

#%%
models = ["MPI_GE_onepct", "MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]

all_model_data = {
    model: read_eof_decade(model) for model in models
}
#%%
# 90th and 10th quantile of all models
quantile_90 = {
    model: quantile_90_10(all_model_data[model], 0.9) for model in models
}

quantile_10 = {
    model: quantile_90_10(all_model_data[model], 0.1) for model in models
}

#%%
# Percentage of values <= 1.5
percentage_1_5 = {
    model: percentage_le_1_5(all_model_data[model], 1.5) for model in models
}

percentage_minus_1_5 = {
    model: percentage_le_1_5(all_model_data[model], -1.5) for model in models
}
# %%

# plotting
fig, ax = plt.subplots(1, 2, figsize=(12, 6))

# Plotting 90th and 10th quantiles
bars_90 = ax[0].bar(models, [quantile_90[model] for model in models], label="90th quantile")
bars_10 = ax[0].bar(models, [quantile_10[model] for model in models], label="10th quantile", alpha=0.7)
ax[0].set_xlabel("Models")
ax[0].set_ylabel("Quantiles")
ax[0].legend()
ax[0].set_ylim(-2, 2)

# Adding value labels on top of bars
for bar in bars_90:
    yval = bar.get_height()
    ax[0].text(bar.get_x() + bar.get_width()/2, yval, f'{yval:.1f}', ha='center', va='bottom')

for bar in bars_10:
    yval = bar.get_height()
    ax[0].text(bar.get_x() + bar.get_width()/2, yval - 0.1, f'{yval:.1f}', ha='center', va='bottom')

# Plotting percentage of values <= 1.5 and -1.5
bars_1_5 = ax[1].bar(models, [percentage_1_5[model] for model in models], label="<= 1.5")
bars_minus_1_5 = ax[1].bar(models, [percentage_minus_1_5[model] for model in models], label="<= -1.5", alpha=0.7)
ax[1].set_xlabel("Models")
ax[1].set_ylabel("Percentage")
ax[1].legend()
ax[1].set_ylim(0, 110)

# Adding value labels on top of bars
for bar in bars_1_5:
    yval = bar.get_height()
    ax[1].text(bar.get_x() + bar.get_width()/2, yval, f'{yval:.1f}', ha='center', va='bottom')

for bar in bars_minus_1_5:
    yval = bar.get_height()
    ax[1].text(bar.get_x() + bar.get_width()/2, yval, f'{yval:.1f}', ha='center', va='bottom')

plt.tight_layout()
plt.show()


# %%
quantile_93 = {
    model: quantile_90_10(all_model_data[model], 0.93) for model in models
}

quantile_7 = {
    model: quantile_90_10(all_model_data[model], 0.07) for model in models
}
# %%

# plotting
fig, ax = plt.subplots(1, 2, figsize=(12, 6))

# Plotting 90th and 10th quantiles
bars_90 = ax[0].bar(models, [quantile_90[model] for model in models], label="90th quantile")
bars_10 = ax[0].bar(models, [quantile_10[model] for model in models], label="10th quantile", alpha=0.7)
ax[0].set_xlabel("Models")
ax[0].set_ylabel("Quantiles")
ax[0].legend()
ax[0].set_ylim(-2, 2)

# Adding value labels on top of bars
for bar in bars_90:
    yval = bar.get_height()
    ax[0].text(bar.get_x() + bar.get_width()/2, yval, f'{yval:.1f}', ha='center', va='bottom')

for bar in bars_10:
    yval = bar.get_height()
    ax[0].text(bar.get_x() + bar.get_width()/2, yval - 0.1, f'{yval:.1f}', ha='center', va='bottom')

# Plotting percentage of values <= 1.5 and -1.5
bars_1_5 = ax[1].bar(models, [percentage_1_5[model] for model in models], label="<= 1.5")
bars_minus_1_5 = ax[1].bar(models, [percentage_minus_1_5[model] for model in models], label="<= -1.5", alpha=0.7)
ax[1].set_xlabel("Models")
ax[1].set_ylabel("Percentage")
ax[1].legend()
ax[1].set_ylim(0, 110)

# Adding value labels on top of bars
for bar in bars_1_5:
    yval = bar.get_height()
    ax[1].text(bar.get_x() + bar.get_width()/2, yval, f'{yval:.1f}', ha='center', va='bottom')

for bar in bars_minus_1_5:
    yval = bar.get_height()
    ax[1].text(bar.get_x() + bar.get_width()/2, yval, f'{yval:.1f}', ha='center', va='bottom')

# Plotting 93rd and 7th quantiles as lines
ax[0].plot(models, [quantile_93[model] for model in models], label="93rd quantile", color='red', marker='o')
ax[0].plot(models, [quantile_7[model] for model in models], label="7th quantile", color='blue', marker='o')
ax[0].legend(ncol = 2)

# add a b
ax[0].text(0.01, 0.98, 'a', transform=ax[0].transAxes, fontsize=16, fontweight='bold', va='top', ha='left')
ax[1].text(0.01, 0.98, 'b', transform=ax[1].transAxes, fontsize=16, fontweight='bold', va='top', ha='left')


plt.tight_layout()
plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/paper_supplymentary/percentage_std.png")

# %%
