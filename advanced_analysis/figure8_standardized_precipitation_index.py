import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from matplotlib.ticker import MultipleLocator

# === File paths ===
folder_path = '/data/merra2_prev_nwc'
file_paths = [os.path.join(folder_path, f'prev_{year}.nc') for year in range(1982, 2021)]

# === Calculate weighted average precipitation for each year ===
pre_int_mean = []
for file_path in file_paths:
    ds = xr.open_dataset(file_path)
    data = ds['tp'] * 1000  # Convert to mm
    weights = np.cos(np.deg2rad(data.lat))  # Latitude weighting
    weighted_avg = data.weighted(weights).mean(dim=["lat", "lon"])
    pre_int_mean.append(float(weighted_avg))

# === Calculate SPI ===
pre_array = np.array(pre_int_mean)
mu = np.mean(pre_array)
sigma = np.std(pre_array, ddof=1)  # Sample standard deviation
spi = (pre_array - mu) / sigma

# === Output SPI values ===
print("=== SPI Value List ===")
for year, val in zip(range(1982, 2021), spi):
    print(f"{year} SPI: {val:.2f}")

# === Plot bar chart ===
colors = ['red' if val <= -1 else 'blue' if val >= 1 else 'gray' for val in spi]
years = np.arange(1982, 2021)

plt.figure(figsize=(12, 8))
plt.bar(years, spi, color=colors)
plt.axhline(0, color='black', linestyle='--')
plt.axhline(1, color='blue', linestyle='--', linewidth=0.5)
plt.axhline(-1, color='red', linestyle='--', linewidth=0.5)
plt.xlabel('Year', fontsize=16)
plt.ylabel('Standardized Precipitation Variability',fontsize=16)
plt.gca().xaxis.set_minor_locator(MultipleLocator(1))
plt.gca().xaxis.set_major_locator(MultipleLocator(5))
#plt.title('Standardized Precipitation Variability', loc='left')
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.xticks(fontsize=14)  # Adjust x-axis tick label font size
plt.yticks(fontsize=14)  # Adjust y-axis tick label font size

# === Add red annotation for dry years (SPI <= -1) in bottom right ===
dry_years = [str(year) for year, val in zip(years, spi) if val <= -1]
if dry_years:
    dry_text = 'Low-variability years: ' + ', '.join(dry_years)
    plt.text(0.98, 0.02, dry_text, transform=plt.gca().transAxes,
             color='red', ha='right', va='bottom', fontsize=20)

# === Add blue annotation for wet years (SPI >= 1) in top left ===
wet_years = [str(year) for year, val in zip(years, spi) if val >= 1]
if wet_years:
    wet_text = 'High-variability years: ' + ', '.join(wet_years)
    plt.text(0.02, 0.98, wet_text, transform=plt.gca().transAxes,
             color='blue', ha='left', va='top', fontsize=20)

# === Save figure ===
plt.savefig('', dpi=800)
plt.show()