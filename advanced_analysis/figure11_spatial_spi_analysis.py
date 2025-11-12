import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmaps
import cartopy.io.shapereader as shpreader
import numpy as np
from matplotlib.patches import Rectangle
import matplotlib.ticker as ticker
from matplotlib.colors import BoundaryNorm

fig = plt.figure(figsize=(28, 26))
def format_func(value, tick_number):
    return f'{value * 1000:.1f}'

# ------------------------ (a)----------------------
file_path = '/data1/tangzhen/moisture_flux_mean/wet_year_mositure_flux.nc'  # Dry years or wet years
title = 'High-variability years'

# === Read data ===
ds = xr.open_dataset(file_path)
ds_subset = ds.sel(lon=slice(50, 160), lat=slice(0, 80))
uq = ds_subset['qu']
vq = ds_subset['qv']
div = ds_subset['divQ']*10000

# === Sparse grid ===
lon_sparse = ds_subset.lon[::1]
lat_sparse = ds_subset.lat[::1]
uq_sparse = uq[::1, ::1]
vq_sparse = vq[::1, ::1]

# === Calculate flux intensity magnitude (for arrow legend) ===
magnitude = np.sqrt(uq**2 + vq**2)
mean_flux = float(magnitude.mean())
print(f'Average water vapor flux intensity: {mean_flux:.2f} kg/m/s')

# === Set arrow legend value (uniformly set to 120) ===
arrow_value = 135

# === Plotting ===
ax_a = fig.add_subplot(3, 1, 1, projection=ccrs.PlateCarree())

# Read and draw shp
china = shpreader.Reader('/home/tangzhen/northwest/northwest_China.shp').geometries()
ax_a.add_geometries(china, ccrs.PlateCarree(), facecolor='none', edgecolor='green', linewidth=3, zorder=2)

# Add geographic features
ax_a.coastlines()
ax_a.add_feature(cfeature.BORDERS, linestyle=':')
ax_a.add_feature(cfeature.LAND, facecolor='lightgray')
ax_a.add_feature(cfeature.OCEAN, facecolor='white')
ax_a.set_extent([50, 150, 00, 60], crs=ccrs.PlateCarree())

# Set longitude and latitude ticks
xticks = np.arange(60, 160, 10)  # Remove 60 and 120
yticks = np.arange(10, 80, 10)   # Remove 20 and 60
ax_a.set_xticks(xticks, crs=ccrs.PlateCarree())
ax_a.set_yticks(yticks, crs=ccrs.PlateCarree())
ax_a.set_xticks(xticks, crs=ccrs.PlateCarree())
ax_a.set_yticks(yticks, crs=ccrs.PlateCarree())

# Set tick label format
ax_a.set_xticklabels([f'{lon}°E' for lon in xticks], fontsize=18)
ax_a.set_yticklabels([f'{lat}°N' for lat in yticks], fontsize=18)

# Show tick lines (ticks)
ax_a.tick_params(axis='both', which='major', length=5, width=1, labelsize=10)


# === Calculate maximum and minimum values of data ===
div_max = np.max(div.values)  # Calculate maximum value
div_min = np.min(div.values)  # Calculate minimum value

# Output maximum and minimum values
print(f"div_max: {div_max}, div_min: {div_min}")
boundaries1 = np.arange(-2.5, 0, 0.25)
boundaries2 = np.arange(0, 2.51, 0.25)
boundaries = np.concatenate([boundaries1, boundaries2])
norm = BoundaryNorm(boundaries, ncolors=64)

# === Draw divergence background ===
div_plot = ax_a.contourf(
    ds_subset.lon, ds_subset.lat, div,
    transform=ccrs.PlateCarree(),
    cmap=cmaps.cmp_b2r,
    levels=boundaries,  # Use symmetric color bar range
    norm = norm,
    extend='both',
)

# === Draw uniformly scaled vector arrows ===
q = ax_a.quiver(
    lon_sparse, lat_sparse, uq_sparse, vq_sparse,
    transform=ccrs.PlateCarree(),
    width=0.003,
    regrid_shape=20,
    scale=3500 # Uniform scaling ratio
)

# === Add arrow reference legend ===
ax_a.quiverkey(
    q, 1.1, 0.05, arrow_value,
    f'{arrow_value} kg·m⁻²·s⁻¹',
    labelpos='E', coordinates='axes',
    fontproperties={'size': 12}
)

# === Add average flux annotation ===
# === Add title and color bar ===
ax_a.set_title(title, fontsize=20)
cb_a = plt.colorbar(div_plot, ax=ax_a, orientation='vertical',
             shrink=0.7,      # Control length (between 0-1, smaller is shorter)
             aspect=20)       # Control width (default is 20, smaller is thinner)

# If you need larger tick fonts, you can use
cb_a.ax.tick_params(labelsize=12)
cb_a.ax.annotate(r'$\times 10^{-4}$', xy=(0.5, 1.08), xycoords='axes fraction',
                 ha='center', va='center', fontsize=12, fontweight='bold')
cb_a.set_label('kg·m⁻²·s⁻¹', fontsize=16)
# === Add red rectangle between longitude/latitude range (105 to 112, 30 to 40) ===
rect = Rectangle((105, 30), 6, 8, linewidth=3, edgecolor='red', facecolor='none', linestyle='-')
ax_a.add_patch(rect)
ax_a.text(0.00, 1.01, "(a)", transform=plt.gca().transAxes, fontsize=18, ha='left', va='bottom', fontweight='bold')

# ------------------------ (b)----------------------
file_path = '/data1/tangzhen/moisture_flux_mean/dry_year_mositure_flux.nc'  # Dry years or wet years
title = 'Low-variability years'

# === Read data ===
ds = xr.open_dataset(file_path)
ds_subset = ds.sel(lon=slice(50, 160), lat=slice(0, 80))
uq = ds_subset['qu']
vq = ds_subset['qv']
div = ds_subset['divQ']*10000

# === Sparse grid ===
lon_sparse = ds_subset.lon[::1]
lat_sparse = ds_subset.lat[::1]
uq_sparse = uq[::1, ::1]
vq_sparse = vq[::1, ::1]

# === Calculate flux intensity magnitude (for arrow legend) ===
magnitude = np.sqrt(uq**2 + vq**2)
mean_flux = float(magnitude.mean())
print(f'Average water vapor flux intensity: {mean_flux:.2f} kg/m/s')

# === Set arrow legend value (uniformly set to 120) ===
arrow_value = 135  

# === Plotting ===
ax_b = fig.add_subplot(3, 1, 2, projection=ccrs.PlateCarree())

# Read and draw shp
china = shpreader.Reader('/home/tangzhen/northwest/northwest_China.shp').geometries()
ax_b.add_geometries(china, ccrs.PlateCarree(), facecolor='none', edgecolor='green', linewidth=3, zorder=2)

# Add geographic features
ax_b.coastlines()
ax_b.add_feature(cfeature.BORDERS, linestyle=':')
ax_b.add_feature(cfeature.LAND, facecolor='lightgray')
ax_b.add_feature(cfeature.OCEAN, facecolor='white')
ax_b.set_extent([50, 150, 00, 60], crs=ccrs.PlateCarree())

# Set longitude and latitude ticks
xticks = np.arange(60, 160, 10)  # Remove 60 and 120
yticks = np.arange(10, 80, 10)   # Remove 20 and 60
ax_b.set_xticks(xticks, crs=ccrs.PlateCarree())
ax_b.set_yticks(yticks, crs=ccrs.PlateCarree())
ax_b.set_xticks(xticks, crs=ccrs.PlateCarree())
ax_b.set_yticks(yticks, crs=ccrs.PlateCarree())

# Set tick label format
ax_b.set_xticklabels([f'{lon}°E' for lon in xticks], fontsize=18)
ax_b.set_yticklabels([f'{lat}°N' for lat in yticks], fontsize=18)

# Show tick lines (ticks)
ax_b.tick_params(axis='both', which='major', length=5, width=1, labelsize=10)


# === Calculate maximum and minimum values of data ===
div_max = np.max(div.values)  # Calculate maximum value
div_min = np.min(div.values)  # Calculate minimum value

# Output maximum and minimum values
print(f"div_max: {div_max}, div_min: {div_min}")
boundaries1 = np.arange(-2.5, 0, 0.25)
boundaries2 = np.arange(0, 2.51, 0.25)
boundaries = np.concatenate([boundaries1, boundaries2])
norm = BoundaryNorm(boundaries, ncolors=64)

# === Draw divergence background ===
div_plot = ax_b.contourf(
    ds_subset.lon, ds_subset.lat, div,
    transform=ccrs.PlateCarree(),
    cmap=cmaps.cmp_b2r,
    levels=boundaries, 
    norm = norm, # Use symmetric color bar range
    extend='both',
)

# === Draw uniformly scaled vector arrows ===
q = ax_b.quiver(
    lon_sparse, lat_sparse, uq_sparse, vq_sparse,
    transform=ccrs.PlateCarree(),
    width=0.003,
    regrid_shape=20,
    scale=3500 # Uniform scaling ratio
)

# === Add arrow reference legend ===
ax_b.quiverkey(
    q, 1.1, 0.05, arrow_value,
    f'{arrow_value} kg·m⁻²·s⁻¹',
    labelpos='E', coordinates='axes',
    fontproperties={'size': 12}
)


ax_b.set_title(title, fontsize=20)
cb_b = plt.colorbar(div_plot, ax=ax_b, orientation='vertical',
             shrink=0.7,      # Control length (between 0-1, smaller is shorter)
             aspect=20)       # Control width (default is 20, smaller is thinner)

# If you need larger tick fonts, you can use
cb_b.ax.tick_params(labelsize=12)
cb_b.ax.annotate(r'$\times 10^{-4}$', xy=(0.5, 1.08), xycoords='axes fraction',
                 ha='center', va='center', fontsize=12, fontweight='bold')
cb_b.set_label('kg·m⁻²·s⁻¹', fontsize=16)
# === Add red rectangle between longitude/latitude range (105 to 112, 30 to 40) ===
rect = Rectangle((105, 30), 6, 8, linewidth=3, edgecolor='red', facecolor='none', linestyle='-')
ax_b.add_patch(rect)
ax_b.text(0.00, 1.01, "(b)", transform=plt.gca().transAxes, fontsize=18, ha='left', va='bottom', fontweight='bold')

# ------------------------ (c)----------------------
# === Modify here ===
file_path = '/data1/tangzhen/moisture_flux_mean/difference_mositure_flux.nc'  # Dry years or wet years
title = 'Difference'

# === Read data ===
ds = xr.open_dataset(file_path)
ds_subset = ds.sel(lon=slice(50, 160), lat=slice(0, 80))
uq = ds_subset['qu']
vq = ds_subset['qv']
div = ds_subset['divQ']*10000

# === Sparse grid ===
lon_sparse = ds_subset.lon[::1]
lat_sparse = ds_subset.lat[::1]
uq_sparse = uq[::1, ::1]
vq_sparse = vq[::1, ::1]

# === Calculate flux intensity magnitude (for arrow legend) ===
magnitude = np.sqrt(uq**2 + vq**2)
mean_flux = float(magnitude.mean())
print(f'Average water vapor flux intensity: {mean_flux:.2f} kg/m/s')

# === Set arrow legend value (uniformly set to 120) ===
arrow_value = 35  # Optional 100, 120, 130, uniform for dry and wet years for comparison

# === Plotting ===
ax = fig.add_subplot(3, 1, 3, projection=ccrs.PlateCarree())


# Read and draw shp
china = shpreader.Reader('/home/tangzhen/northwest/northwest_China.shp').geometries()
ax.add_geometries(china, ccrs.PlateCarree(), facecolor='none', edgecolor='green', linewidth=3, zorder=2)

# Add geographic features
ax.coastlines()
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')
ax.set_extent([50, 150, 00, 60], crs=ccrs.PlateCarree())

# Set longitude and latitude ticks
xticks = np.arange(60, 160, 10)  # Remove 60 and 120
yticks = np.arange(10, 80, 10)   # Remove 20 and 60
ax.set_xticks(xticks, crs=ccrs.PlateCarree())
ax.set_yticks(yticks, crs=ccrs.PlateCarree())
ax.set_xticks(xticks, crs=ccrs.PlateCarree())
ax.set_yticks(yticks, crs=ccrs.PlateCarree())

# Set tick label format
ax.set_xticklabels([f'{lon}°E' for lon in xticks], fontsize=18)
ax.set_yticklabels([f'{lat}°N' for lat in yticks], fontsize=18)

# Show tick lines (ticks)
ax.tick_params(axis='both', which='major', length=5, width=1, labelsize=10)


# === Calculate maximum and minimum values of data ===
div_max = np.max(div.values)  # Calculate maximum value
div_min = np.min(div.values)  # Calculate minimum value

# Output maximum and minimum values
print(f"div_max: {div_max}, div_min: {div_min}")
boundaries1 = np.arange(-1, 0, 0.1)
boundaries2 = np.arange(0, 1.1, 0.1)
boundaries = np.concatenate([boundaries1, boundaries2])
norm = BoundaryNorm(boundaries, ncolors=64)

# === Draw divergence background ===
div_plot = ax.contourf(
    ds_subset.lon, ds_subset.lat, div,
    transform=ccrs.PlateCarree(),
    cmap=cmaps.cmp_b2r,
    levels=boundaries, 
    norm = norm, # Use symmetric color bar range
    extend='both',
)
# === Draw uniformly scaled vector arrows ===
q = ax.quiver(
    lon_sparse, lat_sparse, uq_sparse, vq_sparse,
    transform=ccrs.PlateCarree(),
    width=0.003,
    regrid_shape=20,
    scale=800# Uniform scaling ratio
)


ax.quiverkey(
    q, 1.1, 0.05, arrow_value,
    f'{arrow_value} kg·m⁻²·s⁻¹',
    labelpos='E', coordinates='axes',
    fontproperties={'size': 12}
)

# === Add average flux annotation ===
# === Add title and color bar ===
ax.set_title(title, fontsize=20)
cbar = plt.colorbar(div_plot, ax=ax, orientation='vertical',
             shrink=0.7,      # Control length (between 0-1, smaller is shorter)
             aspect=20)       # Control width (default is 20, smaller is thinner)

# If you need larger tick fonts, you can use
cbar.ax.tick_params(labelsize=12)
cbar.ax.annotate(r'$\times 10^{-4}$', xy=(0.5, 1.08), xycoords='axes fraction',
                 ha='center', va='center', fontsize=12, fontweight='bold')
cbar.set_label('kg·m⁻²·s⁻¹', fontsize=16)
rect = Rectangle((105, 30), 6, 8, linewidth=3, edgecolor='red', facecolor='none', linestyle='-')
ax.add_patch(rect)
ax.text(0.00, 1.01, "(c)", transform=plt.gca().transAxes, fontsize=18, ha='left', va='bottom', fontweight='bold')


# === Save figure ===
plt.savefig('/home/tangzhen/Result/moisture_flux_np.png', dpi=600, bbox_inches="tight")