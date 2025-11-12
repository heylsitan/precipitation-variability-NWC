import os
import pandas as pd
import numpy as np
import geopandas as gpd
import xarray as xr
from scipy.stats import linregress, pearsonr
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy.spatial import Voronoi
from shapely.geometry import Point, Polygon
from sklearn.neighbors import NearestNeighbors

# ====================== Station Data Processing ======================

station_df = pd.read_excel('/data/station_prev.xlsx')
station_years = station_df['Year'].values
station_mean = station_df['Prev'].values
station_slope, station_intercept, station_r, station_p, station_err = linregress(station_years, station_mean)
station_trend = station_slope * station_years + station_intercept

# Keep only 1982–2018 to align with MERRA‑2
station_years = station_years[station_years <= 2018]
station_mean = station_mean[:len(station_years)]  # Ensure lengths match

# ====================== MERRA‑2 Data Processing ======================

merra_folder = '/data//prev/merra2_prev_nwc'
# Keep 1982–2018 for correlation calculation
merra_years_for_corr = np.arange(1982, 2019)  # 1982-2018
merra_files_for_corr = [os.path.join(merra_folder, f'prev_{year}.nc') for year in merra_years_for_corr]
merra_mean_for_corr = []

for file in merra_files_for_corr:
    ds = xr.open_dataset(file)
    data = ds['tp'] * 1000  # Convert to mm
    weights = np.cos(np.deg2rad(data.lat))
    weighted_avg = data.weighted(weights).mean(dim=["lat", "lon"])
    merra_mean_for_corr.append(float(weighted_avg))

# Keep only 1982–2018
merra_years_for_corr = merra_years_for_corr[merra_years_for_corr <= 2018]
merra_mean_for_corr = merra_mean_for_corr[:len(merra_years_for_corr)]  # Ensure lengths match

merra_slope_corr, merra_intercept_corr, merra_r_corr, merra_p_full, merra_err_full = linregress(merra_years_for_corr, merra_mean_for_corr)

# Compute Pearson correlation
correlation, p_value = pearsonr(station_mean, merra_mean_for_corr)
print(f"Pearson correlation (1982-2018): {correlation}, p-value: {p_value}")

# ====================== MERRA‑2 Processing — Full Data for Trend Plotting ======================

# Use 1982–2020 for plotting and trend computation
merra_years_full = np.arange(1982, 2021)  # 1982-2020
merra_files_full = [os.path.join(merra_folder, f'prev_{year}.nc') for year in merra_years_full]
merra_mean_full = []

for file in merra_files_full:
    ds = xr.open_dataset(file)
    data = ds['tp'] * 1000  # Convert to mm
    weights = np.cos(np.deg2rad(data.lat))
    weighted_avg = data.weighted(weights).mean(dim=["lat", "lon"])
    merra_mean_full.append(float(weighted_avg))

# Perform linear regression
merra_slope_full, merra_intercept_full, merra_r_full, merra_p_full, merra_err_full = linregress(merra_years_full, merra_mean_full)
merra_trend_full = merra_slope_full * merra_years_full + merra_intercept_full

# ====================== Visualization ======================

fig, axs = plt.subplots(1, 2, figsize=(16, 6), sharey=False)
plt.subplots_adjust(wspace=0.50)  # wspace is horizontal spacing; default 0.2

# --- Left panel: Station ---
axs[0].plot(station_years, station_mean, marker='s', color='r', markersize=8, linewidth=2, linestyle='--')
axs[0].plot(station_years, station_trend, linestyle='-', color='black', linewidth=2)
axs[0].annotate(f'Trend={station_slope:.4f}mm/year(p<0.05)',
                xy=(station_years[0], station_trend[0]),
                xytext=(-10, 215),
                textcoords='offset points',
                color='r',
                fontsize=16)
axs[0].text(0.01, 1.01, "(a)", transform=axs[0].transAxes, fontsize=18, fontweight='bold')
axs[0].set_title('Station', fontsize=18)
axs[0].set_xlabel('Years', fontsize=16)
axs[0].set_ylabel('mm', fontsize=16)

# --- Right panel: MERRA‑2 ---
axs[1].plot(merra_years_full, merra_mean_full, marker='s', color='blue', markersize=8, linewidth=2, linestyle='--')
axs[1].plot(merra_years_full, merra_trend_full, linestyle='-', color='black', linewidth=2)

# Add correlation text to right panel (b)
axs[1].annotate(f'Trend(1982-2020)={merra_slope_full:.4f}mm/year(p<0.05)\n'
                f'Trend(1982-2018)={merra_slope_corr:.4f}mm/year(p<0.05)\n'
                f'R(1982-2018)={correlation:.2f}',
                xy=(merra_years_full[0], merra_trend_full[0]),
                xytext=(-10, 222),
                textcoords='offset points',
                color='blue',
                fontsize=16)
axs[1].text(0.01, 1.01, "(b)", transform=axs[1].transAxes, fontsize=18, fontweight='bold')
axs[1].set_title('MERRA2', fontsize=18)
axs[1].set_xlabel('Years', fontsize=16)
axs[1].set_ylabel('mm', fontsize=16)

# Axis styling
for ax in axs:
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_major_locator(MultipleLocator(0.2))
    ax.tick_params(axis='both', which='minor', width=1, length=3, labelsize=16)
    ax.tick_params(axis='both', which='major', width=1, length=6, labelsize=16)

plt.tight_layout()
plt.savefig('', dpi=800, bbox_inches='tight')
plt.savefig('', bbox_inches='tight')
plt.show()
