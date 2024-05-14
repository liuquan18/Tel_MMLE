#%%
# print hello world
print("Hello World!")

#%%
print("hello my copolit")
# that's it!import numpy as np

#%%
import statsmodels.api as sm
import numpy as np
import matplotlib.pyplot as plt
#%%
# Set the parameters of the AR(1) model
phi = 0.3
sigma = 1.0
n = 1000

# Generate a time series from the AR(1) model
np.random.seed(123)
epsilon = np.random.normal(loc=0, scale=sigma, size=n)
y = np.zeros(n)
y[0] = epsilon[0]
for i in range(1, n):
    y[i] = phi * y[i-1] + epsilon[i]

# Fit an AR(1) model to the time series
model = sm.tsa.ARIMA(y, order=(1, 0, 0))
results = model.fit()

# Print the model summary
print(results.summary())
# %%
# Generate a series from the AR(1) model
series = results.predict(start=0, end=n-1)

# Print the series
print(series)
# %%
series.mean()
# %%
series.std()
# %%
plt.hist(series, bins=20)
# %%
plt.plot(series)

#%%
# Set the parameters of the AR(1) model
phi = 0.3
sigma = 1.0
n = 1000
n_realizations = 9000

# Generate the realizations
np.random.seed(123)
epsilon = np.random.normal(loc=0, scale=sigma, size=(n_realizations, n))
y = np.zeros((n_realizations, n))
y[:, 0] = epsilon[:, 0]
for i in range(1, n):
    y[:, i] = phi * y[:, i-1] + epsilon[:, i]

# Compute the mean and standard deviation of the realizations
mean = np.mean(y, axis=0)
std = np.std(y, axis=0)

# Compute the 5%-95% confidence interval
lower = np.percentile(y, 5, axis=0)
upper = np.percentile(y, 95, axis=0)

# Plot the mean and confidence interval
plt.plot(mean, color='blue', label='Mean')
plt.fill_between(range(n), lower, upper, color='gray', alpha=0.5, label='5%-95% CI')
plt.legend()
plt.show()# %%

# %%
# Count the number of records above 1.5 and below -1.5
above = np.sum(y[0] > 1.5)
below = np.sum(y[0] < -1.5)

# Print the results
print(f"Number of records above 1.5: {above}")
print(f"Number of records below -1.5: {below}")
# %%
# Count the number of records above 1.5 and below -1.5 in all realizations
above_all = np.sum(y > 1.5, axis=1)
below_all = np.sum(y < -1.5, axis=1)

# Compute the 5%-95% confidence interval of the counts
above_lower, above_upper = np.percentile(above_all, [5, 95])
below_lower, below_upper = np.percentile(below_all, [5, 95])

# Print the results
print(f"5%-95% interval of records above 1.5: ({above_lower}, {above_upper})")
print(f"5%-95% interval of records below -1.5: ({below_lower}, {below_upper})")
# %%
# Set the parameters of the AR(1) model
sigma = 1.0
n = 1000
n_realizations = 9000
phis = np.linspace(0, 1, 21)

# Compute the confidence intervals for different values of phi
intervals = []
for phi in phis:
    # Generate the realizations
    np.random.seed(123)
    epsilon = np.random.normal(loc=0, scale=sigma, size=(n_realizations, n))
    y = np.zeros((n_realizations, n))
    y[:, 0] = epsilon[:, 0]
    for i in range(1, n):
        y[:, i] = phi * y[:, i-1] + epsilon[:, i]

    # Compute the 5%-95% confidence interval
    lower = np.percentile(y, 5, axis=0)
    upper = np.percentile(y, 95, axis=0)
    intervals.append((lower, upper))

# Plot the confidence intervals
fig, ax = plt.subplots()
for i, phi in enumerate(phis):
    ax.plot(intervals[i][0], color='gray', alpha=0.5)
    ax.plot(intervals[i][1], color='gray', alpha=0.5)
    ax.fill_between(range(n), intervals[i][0], intervals[i][1], alpha=0.2, label=f'phi={phi:.1f}')
ax.legend()
plt.show()
# %%
# Set the parameters of the AR(1) model
sigma = 1.0
n = 1000
n_realizations = 9000
phis = np.linspace(0, 1, 101)

# Compute the upper boundary of the confidence interval for records above 1.5
upper_bounds = []
for phi in phis:
    # Generate the realizations
    np.random.seed(123)
    epsilon = np.random.normal(loc=0, scale=sigma, size=(n_realizations, n))
    y = np.zeros((n_realizations, n))
    y[:, 0] = epsilon[:, 0]
    
    for i in range(1, n):
        y[:, i] = phi * y[:, i-1] + epsilon[:, i]

    # Compute the 95% confidence interval for records above 1.5
    upper = np.percentile(y[y > 1.5], 95,axis=0)
    upper_bounds.append(upper)
#%%
# Plot the upper boundary of the confidence interval
plt.plot(phis, upper_bounds)
plt.xlabel('phi')
plt.ylabel('Upper boundary of 95% CI for records above 1.5')
plt.xlim(0,0.6)
plt.show()
# %%
