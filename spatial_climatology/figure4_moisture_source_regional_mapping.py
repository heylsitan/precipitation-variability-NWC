import numpy as np
import cmaps
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from netCDF4 import Dataset
from matplotlib.colors import BoundaryNorm
from matplotlib import colors

fig = plt.figure(figsize=(16, 12))  # Horizontal layout, suitable for two figures + colorbar

# === (a) e_track_contourf ===
file_path = '/data/wam2_result/mean_1982_2020.nc'
nc_data = Dataset(file_path, 'r')
var_name = 'e_track'
variable = nc_data.variables[var_name][:]
lon = nc_data.variables['longitude'][:]
lat = nc_data.variables['latitude'][:]

# Main figure a: Map
ax_a = fig.add_axes([0.06, 0.1, 0.38, 0.8], projection=ccrs.PlateCarree())
ax_a.set_extent([-30, 180, -50, 90], crs=ccrs.PlateCarree())
ax_a.coastlines()

xticks = np.arange(0, 180, 30)
yticks = np.arange(-40, 90, 20)
ax_a.set_xticks(xticks)
ax_a.set_yticks(yticks)
ax_a.set_xticklabels([f'{x}°E' for x in xticks], fontsize=16)
yticklabels = [f'{abs(y)}°{"N" if y >= 0 else "S"}' for y in yticks]
ax_a.set_yticklabels(yticklabels, fontsize=16)
ax_a.grid(True, color='gray', linestyle='-', linewidth=0.5)

# Set color breakpoints
boundariesf = np.arange(0, 0.1)
boundaries0 = np.arange(0.1, 0.5, 0.2)  # Finer intervals
boundaries1 = np.arange(0.5, 8.0, 1.5)
boundaries2 = np.arange(8.0, 24.0 ,4.0)
boundaries3 = np.arange(24.0, 80.0, 12.0)
boundaries = np.concatenate((boundariesf, boundaries0, boundaries1, boundaries2, boundaries3))
norm = BoundaryNorm(boundaries, ncolors=18)

tp_plot = ax_a.contourf(lon, lat, variable, cmap=cmaps.precip3_16lev,
                        extend='max', levels=boundaries, norm=norm, transform=ccrs.PlateCarree())

# Add China boundary
china_shp = '/home/tangzhen/northwest/northwest_China.shp'
china_geom = shpreader.Reader(china_shp).geometries()
ax_a.add_geometries(china_geom, ccrs.PlateCarree(),
                    facecolor='none', edgecolor='red', linewidth=1.2, zorder=3)

# Add colorbar, without compressing the main figure
cax = fig.add_axes([0.45, 0.33, 0.015, 0.33])  # Manually add colorbar to the right of figure a
cb_a = plt.colorbar(tp_plot, cax=cax, orientation='vertical')
cb_a.set_ticks(boundaries[:-1])
cb_a.ax.tick_params(which='minor', length=0)
cb_a.set_label('mm', fontsize=16)
cb_a.ax.tick_params(labelsize=12)

ax_a.set_title("Moisture source", fontsize=18, loc='center')
ax_a.text(0.00, 1.02, "(a)", transform=ax_a.transAxes, fontsize=18, ha='left', va='bottom', fontweight='bold')


# === (b) source_region ===
ax_b = fig.add_axes([0.56, 0.1, 0.38, 0.8], projection=ccrs.PlateCarree())
ax_b.set_extent([-30, 180, -50, 90], crs=ccrs.PlateCarree())
ax_b.coastlines()
ax_b.add_feature(cfeature.LAND, facecolor='lightgray')

ax_b.set_xticks(xticks)
ax_b.set_yticks(yticks)
ax_b.set_xticklabels([f'{x}°E' for x in xticks], fontsize=16)
ax_b.set_yticklabels(yticklabels, fontsize=16)
ax_b.tick_params(axis='x', labelsize=16)
ax_b.tick_params(axis='y', labelsize=16)
ax_b.grid(True, color='gray', linestyle='-', linewidth=0.5)

# Function to add regional layers
def add_region(ax, shp_path, colorcode, label, label_coords):
    geom = shpreader.Reader(shp_path).geometries()
    ax.add_geometries(geom, ccrs.PlateCarree(),
                      facecolor=colors.to_rgba(colorcode, alpha=0.8),
                      edgecolor=colors.to_rgba(colorcode, alpha=0.8),
                      linewidth=1.5, zorder=3)
    ax.text(*label_coords, label, transform=ccrs.PlateCarree(),
            fontsize=18, fontweight='bold', color='black', ha='center')

# Various regions
add_region(ax_b, '/data/shp/NW_region.shp', '#CA0E12', 'NW', (50, 60))
add_region(ax_b, '/data/shp/NE_region.shp', '#25377F', 'NE', (115, 60))
add_region(ax_b, '/data/shp/SW_region.shp', '#2AA7DE', 'SW', (62, 0))
add_region(ax_b, '/data/shp/SE_region.shp', '#F6BD21', 'SE', (130, 15))
add_region(ax_b, '/homenorthwest/northwest_China.shp', 'lightgreen', 'NWC', (90, 38))

ax_b.set_title("Source regions", fontsize=18, loc='center')
ax_b.text(0.00, 1.02, "(b)", transform=ax_b.transAxes, fontsize=18, ha='left', va='bottom', fontweight='bold')

# === Save image ===
plt.savefig('', dpi=600, bbox_inches="tight")