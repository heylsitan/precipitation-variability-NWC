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

station_df = pd.read_excel('')
station_years = station_df['Year'].values
station_mean = station_df['Prev'].values
station_slope, station_intercept, station_r, station_p, station_err = linregress(station_years, station_mean)
station_trend = station_slope * station_years + station_intercept

# Keep only data from 1982 to 2018 to align with MERRA2 data
station_years = station_years[station_years <= 2018]
station_mean = station_mean[:len(station_years)]  # Ensure data length matches

# ====================== MERRA2 Data Processing ======================

merra_folder = ''
merra_years = np.arange(1982, 2021)
merra_files = [os.path.join(merra_folder, f'prev_{year}.nc') for year in merra_years]
merra_mean = []

for file in merra_files:
    ds = xr.open_dataset(file)
    data = ds['tp'] * 1000  # Convert to mm
    weights = np.cos(np.deg2rad(data.lat))
    weighted_avg = data.weighted(weights).mean(dim=["lat", "lon"])
    merra_mean.append(float(weighted_avg))

# Keep only data from 1982 to 2018 to align with station data
merra_years = merra_years[merra_years <= 2018]
merra_mean = merra_mean[:len(merra_years)]  # Ensure data length matches

merra_slope, merra_intercept, merra_r, merra_p, merra_err = linregress(merra_years, merra_mean)
merra_trend = merra_slope * merra_years + merra_intercept

# ====================== Calculate Pearson Correlation Coefficient ======================
correlation, p_value = pearsonr(station_mean, merra_mean)

# ====================== Visualization ======================

fig, axs = plt.subplots(1, 2, figsize=(16, 6), sharey=False)
plt.subplots_adjust(wspace=0.50)  # wspace is the horizontal distance between subplots, default is 0.2

# --- Left plot: Station ---
axs[0].plot(station_years, station_mean, marker='s', color='r', markersize=8, linewidth=2, linestyle='--')
axs[0].plot(station_years, station_trend, linestyle='-', color='black', linewidth=2)
axs[0].annotate(f'Trend={station_slope:.4f}mm/year(p<0.05)',
                xy=(station_years[0], station_trend[0]),
                xytext=(0, 210),
                textcoords='offset points',
                color='r',
                fontsize=18)
axs[0].text(0.01, 1.01, "(a)", transform=axs[0].transAxes, fontsize=18, fontweight='bold')
axs[0].set_title('Station', fontsize=18)
axs[0].set_xlabel('Years', fontsize=16)
axs[0].set_ylabel('mm', fontsize=16)

# --- Right plot: MERRA2 ---
axs[1].plot(merra_years, merra_mean, marker='s', color='blue', markersize=8, linewidth=2, linestyle='--')
axs[1].plot(merra_years, merra_trend, linestyle='-', color='black', linewidth=2)

# Add correlation coefficient to right plot (figure b)
axs[1].annotate(f'Trend={merra_slope:.4f}mm/year(p<0.05)\n'
                f'R={correlation:.2f}',
                xy=(merra_years[0], merra_trend[0]),
                xytext=(0, 230),
                textcoords='offset points',
                color='blue',
                fontsize=18)

axs[1].text(0.01, 1.01, "(b)", transform=axs[1].transAxes, fontsize=18, fontweight='bold')
axs[1].set_title('MERRA2', fontsize=18)
axs[1].set_xlabel('Years', fontsize=16)
axs[1].set_ylabel('mm', fontsize=16)

# Unified coordinate style settings
for ax in axs:
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_major_locator(MultipleLocator(0.2))
    ax.tick_params(axis='both', which='minor', width=1, length=3, labelsize=16)
    ax.tick_params(axis='both', which='major', width=1, length=6, labelsize=16)

plt.tight_layout()
plt.savefig('', dpi=800, bbox_inches='tight')
plt.show()