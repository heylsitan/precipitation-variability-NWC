import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmaps
import cartopy.io.shapereader as shpreader
import numpy as np

fig = plt.figure(figsize=(24, 13))

# ------------------------ (a)----------------------
file_path = '/data1/omage_mean/wet_omage_mean.nc'  # Dry years or wet years
z_file_path = '/data1/z_mean/wet_H_mean.nc'
title = 'High-variability 200 hPa'

# === Read data ===
ds = xr.open_dataset(file_path)
ds_z = xr.open_dataset(z_file_path)
ds_subset = ds.sel(lon=slice(0, 180), lat=slice(0, 80))  # Crop data according to longitude and latitude range
ds_z_subset = ds_z.sel(lon=slice(0, 180), lat=slice(0, 80))

# === Select 200 hPa wind field data ===
u_200hpa = ds_subset['u'].sel(lev=200)  # Select u at 200 hPa level
v_200hpa = ds_subset['v'].sel(lev=200)  # Select v at 200 hPa level
# === Select geopotential height data ===
z_200hpa = ds_z_subset['h'].sel(lev=200)  # Assume geopotential height data is stored in 'h' variable

# === Sparse grid ===
lon_sparse = ds_subset.lon[::3]
lat_sparse = ds_subset.lat[::3]
u_sparse = u_200hpa[::3, ::3]
v_sparse = v_200hpa[::3, ::3]

# === Calculate flux intensity magnitude (for arrow legend) === 
magnitude = np.sqrt(u_200hpa**2 + v_200hpa**2)
mean_flux = float(magnitude.mean())
print(f'Average water vapor flux intensity: {mean_flux:.2f} kg/m/s')

# === Set arrow legend value (uniformly set to 120) === 
arrow_value = 15  # Optional 100, 120, 130, uniform for dry and wet years for comparison

# === Plotting ===
ax_a = fig.add_subplot(2, 2, 1, projection=ccrs.PlateCarree())

# Read and draw shp
china = shpreader.Reader('/home/tangzhen/northwest/northwest_China.shp').geometries()
ax_a.add_geometries(china, ccrs.PlateCarree(), facecolor='none', edgecolor='green', linewidth=1.5, zorder=2, linestyle='-')

# Add geographic features
ax_a.coastlines()
ax_a.add_feature(cfeature.BORDERS, linestyle=':')
ax_a.add_feature(cfeature.LAND, facecolor='none')
ax_a.add_feature(cfeature.OCEAN, facecolor='white')
ax_a.set_extent([60, 120, 20, 60], crs=ccrs.PlateCarree())

# Set longitude and latitude ticks
xticks = np.arange(0, 181, 20)  # Remove 60 and 120
yticks = np.arange(0, 81, 20)   # Remove 20 and 60
ax_a.set_xticks(xticks, crs=ccrs.PlateCarree())
ax_a.set_yticks(yticks, crs=ccrs.PlateCarree())
ax_a.set_xticks(xticks, crs=ccrs.PlateCarree())
ax_a.set_yticks(yticks, crs=ccrs.PlateCarree())

# Set tick label format
ax_a.set_xticklabels([f'{lon}°E' for lon in xticks], fontsize=18)
ax_a.set_yticklabels([f'{lat}°N' for lat in yticks], fontsize=18)

# Show tick lines (ticks)
ax_a.tick_params(axis='both', which='major', length=5, width=1, labelsize=10)

# === Draw uniformly scaled vector arrows ===
q = ax_a.quiver(
    lon_sparse, lat_sparse, u_sparse, v_sparse,
    transform=ccrs.PlateCarree(),
    width=0.002,
    regrid_shape=20,
    scale=400# Uniform scaling ratio
)

# === Add arrow reference legend ===
ax_a.quiverkey(
    q, 0.9, 1.02, arrow_value,
    f'{arrow_value}m/s',
    labelpos='E', coordinates='axes',
    fontproperties={'size': 12}
)

# === Add 12500 geopotential height contour ===
contour_levels = [12500]  # South Asian high corresponds to geopotential height 12500
cs = ax_a.contour(
    ds_z_subset.lon, ds_z_subset.lat, z_200hpa, levels=contour_levels,
    colors='blue', linewidths=2, transform=ccrs.PlateCarree()
)
ax_a.clabel(cs, inline=True, fontsize=10, fmt='%d')

ax_a.text(0.01, 1.01, "(a)", transform=plt.gca().transAxes, fontsize=18, ha='left', va='bottom', fontweight='bold')
# === Add title and color bar ===
ax_a.set_title(title, fontsize=18, loc='center')

# ------------------------ (b)----------------------
file_path = '/data1/omage_mean/dry_omage_mean.nc'  # Dry years or wet years
z_file_path = '/data1/z_mean/dry_H_mean.nc'
title = 'Low-variability 200 hPa'

# === Read data ===
ds = xr.open_dataset(file_path)
ds_z = xr.open_dataset(z_file_path)
ds_subset = ds.sel(lon=slice(0, 180), lat=slice(0, 80))  # Crop data according to longitude and latitude range
ds_z_subset = ds_z.sel(lon=slice(0, 180), lat=slice(0, 80))

# === Select 200 hPa wind field data ===
u_200hpa = ds_subset['u'].sel(lev=200)  # Select u at 200 hPa level
v_200hpa = ds_subset['v'].sel(lev=200)  # Select v at 200 hPa level
# === Select geopotential height data ===
z_200hpa = ds_z_subset['h'].sel(lev=200)  # Assume geopotential height data is stored in 'h' variable

# === Sparse grid ===
lon_sparse = ds_subset.lon[::3]
lat_sparse = ds_subset.lat[::3]
u_sparse = u_200hpa[::3, ::3]
v_sparse = v_200hpa[::3, ::3]

# === Calculate flux intensity magnitude (for arrow legend) === 
magnitude = np.sqrt(u_200hpa**2 + v_200hpa**2)
mean_flux = float(magnitude.mean())
print(f'Average water vapor flux intensity: {mean_flux:.2f} kg/m/s')

# === Set arrow legend value (uniformly set to 120) === 
arrow_value = 15  # Optional 100, 120, 130, uniform for dry and wet years for comparison

# === Plotting ===
ax_b = fig.add_subplot(2, 2, 2, projection=ccrs.PlateCarree())

# Read and draw shp
china = shpreader.Reader('/home/tangzhen/northwest/northwest_China.shp').geometries()
ax_b.add_geometries(china, ccrs.PlateCarree(), facecolor='none', edgecolor='green', linewidth=1.5, zorder=2, linestyle='-')

# Add geographic features
ax_b.coastlines()
ax_b.add_feature(cfeature.BORDERS, linestyle=':')
ax_b.add_feature(cfeature.LAND, facecolor='none')
ax_b.add_feature(cfeature.OCEAN, facecolor='white')
ax_b.set_extent([60, 120, 20, 60], crs=ccrs.PlateCarree())

# Set longitude and latitude ticks
xticks = np.arange(0, 181, 20)  # Remove 60 and 120
yticks = np.arange(0, 81, 20)   # Remove 20 and 60
ax_b.set_xticks(xticks, crs=ccrs.PlateCarree())
ax_b.set_yticks(yticks, crs=ccrs.PlateCarree())
ax_b.set_xticks(xticks, crs=ccrs.PlateCarree())
ax_b.set_yticks(yticks, crs=ccrs.PlateCarree())

# Set tick label format
ax_b.set_xticklabels([f'{lon}°E' for lon in xticks], fontsize=18)
ax_b.set_yticklabels([f'{lat}°N' for lat in yticks], fontsize=18)

# Show tick lines (ticks)
ax_b.tick_params(axis='both', which='major', length=5, width=1, labelsize=10)

# === Draw uniformly scaled vector arrows ===
q = ax_b.quiver(
    lon_sparse, lat_sparse, u_sparse, v_sparse,
    transform=ccrs.PlateCarree(),
    width=0.002,
    regrid_shape=20,
    scale=400# Uniform scaling ratio
)

# === Add arrow reference legend ===
ax_b.quiverkey(
    q, 0.9, 1.02, arrow_value,
    f'{arrow_value}m/s',
    labelpos='E', coordinates='axes',
    fontproperties={'size': 12}
)

# === Add 12500 geopotential height contour ===
contour_levels = [12500]  # South Asian high corresponds to geopotential height 12500
cs = ax_b.contour(
    ds_z_subset.lon, ds_z_subset.lat, z_200hpa, levels=contour_levels,
    colors='blue', linewidths=2, transform=ccrs.PlateCarree()
)
ax_b.clabel(cs, inline=True, fontsize=10, fmt='%d')

ax_b.text(0.01, 1.01, "(b)", transform=plt.gca().transAxes, fontsize=18, ha='left', va='bottom', fontweight='bold')
# === Add title and color bar ===
ax_b.set_title(title, fontsize=18, loc='center')

# ------------------------ (c)----------------------
file_path = '/data1/omage_mean/wet_omage_mean.nc'  # Dry years or wet years
output_path = '/home/Result/wet_WSPH_500hpa.png'  #WSPH
z_file_path = '/data1/z_mean/wet_H_mean.nc'
title = 'High-variability 500 hPa'

# === Read data ===
ds = xr.open_dataset(file_path)
ds_z = xr.open_dataset(z_file_path)
ds_subset = ds.sel(lon=slice(0, 180), lat=slice(0, 80))  # Crop data according to longitude and latitude range
ds_z_subset = ds_z.sel(lon=slice(110, 180), lat=slice(0, 80))

# === Select 200 hPa wind field data ===
u_200hpa = ds_subset['u'].sel(lev=500)  # Select u at 200 hPa level
v_200hpa = ds_subset['v'].sel(lev=500)  # Select v at 200 hPa level
# === Select geopotential height data ===
z_200hpa = ds_z_subset['h'].sel(lev=500)  # Assume geopotential height data is stored in 'h' variable

# === Sparse grid ===
lon_sparse = ds_subset.lon[::3]
lat_sparse = ds_subset.lat[::3]
u_sparse = u_200hpa[::3, ::3]
v_sparse = v_200hpa[::3, ::3]

# === Calculate flux intensity magnitude (for arrow legend) === 
magnitude = np.sqrt(u_200hpa**2 + v_200hpa**2)
mean_flux = float(magnitude.mean())
print(f'Average water vapor flux intensity: {mean_flux:.2f} kg/m/s')

# === Set arrow legend value (uniformly set to 120) === 
arrow_value = 15  # Optional 100, 120, 130, uniform for dry and wet years for comparison

# === Plotting ===
ax_c = fig.add_subplot(2, 2, 3, projection=ccrs.PlateCarree())

# Read and draw shp
china = shpreader.Reader('/home/tangzhen/northwest/northwest_China.shp').geometries()
ax_c.add_geometries(china, ccrs.PlateCarree(), facecolor='none', edgecolor='green', linewidth=1.5, zorder=2, linestyle='-')

# Add geographic features
ax_c.coastlines()
ax_c.add_feature(cfeature.BORDERS, linestyle=':')
ax_c.add_feature(cfeature.LAND, facecolor='none')
ax_c.add_feature(cfeature.OCEAN, facecolor='white')
ax_c.set_extent([60, 120, 20, 60], crs=ccrs.PlateCarree())

# Set longitude and latitude ticks
xticks = np.arange(0, 181, 20)  # Remove 60 and 120
yticks = np.arange(0, 81, 20)   # Remove 20 and 60
ax_c.set_xticks(xticks, crs=ccrs.PlateCarree())
ax_c.set_yticks(yticks, crs=ccrs.PlateCarree())
ax_c.set_xticks(xticks, crs=ccrs.PlateCarree())
ax_c.set_yticks(yticks, crs=ccrs.PlateCarree())

# Set tick label format
ax_c.set_xticklabels([f'{lon}°E' for lon in xticks], fontsize=18)
ax_c.set_yticklabels([f'{lat}°N' for lat in yticks], fontsize=18)

# Show tick lines (ticks)
ax_c.tick_params(axis='both', which='major', length=5, width=1, labelsize=10)

# === Draw uniformly scaled vector arrows ===
q = ax_c.quiver(
    lon_sparse, lat_sparse, u_sparse, v_sparse,
    transform=ccrs.PlateCarree(),
    width=0.002,
    regrid_shape=20,
    scale=400# Uniform scaling ratio
)

# === Add arrow reference legend ===
ax_c.quiverkey(
    q, 0.9, 1.02, arrow_value,
    f'{arrow_value}m/s',
    labelpos='E', coordinates='axes',
    fontproperties={'size': 12}
)

# === Add 12500 geopotential height contour ===
contour_levels = [5880]  # South Asian high corresponds to geopotential height 12500
cs = ax_c.contour(
    ds_z_subset.lon, ds_z_subset.lat, z_200hpa, levels=contour_levels,
    colors='red', linewidths=2, transform=ccrs.PlateCarree()
)
ax_c.clabel(cs, inline=True, fontsize=10, fmt='%d')

# === Add title and color bar ===
ax_c.set_title(title, fontsize=18, loc='center')
ax_c.text(0.01, 1.01, "(c)", transform=plt.gca().transAxes, fontsize=18, ha='left', va='bottom', fontweight='bold')

# ------------------------ (d)----------------------
file_path = '/data1/omage_mean/dry_omage_mean.nc'  # Dry years or wet years
z_file_path = '/data1/z_mean/dry_H_mean.nc'
title = 'Low-variability 500 hPa'

# === Read data ===
ds = xr.open_dataset(file_path)
ds_z = xr.open_dataset(z_file_path)
ds_subset = ds.sel(lon=slice(0, 180), lat=slice(0, 80))  # Crop data according to longitude and latitude range
ds_z_subset = ds_z.sel(lon=slice(110, 180), lat=slice(0, 80))

# === Select 200 hPa wind field data ===
u_200hpa = ds_subset['u'].sel(lev=500)  # Select u at 200 hPa level
v_200hpa = ds_subset['v'].sel(lev=500)  # Select v at 200 hPa level
# === Select geopotential height data ===
z_200hpa = ds_z_subset['h'].sel(lev=500)  # Assume geopotential height data is stored in 'h' variable

# === Sparse grid ===
lon_sparse = ds_subset.lon[::3]
lat_sparse = ds_subset.lat[::3]
u_sparse = u_200hpa[::3, ::3]
v_sparse = v_200hpa[::3, ::3]

# === Calculate flux intensity magnitude (for arrow legend) === 
magnitude = np.sqrt(u_200hpa**2 + v_200hpa**2)
mean_flux = float(magnitude.mean())
print(f'Average water vapor flux intensity: {mean_flux:.2f} kg/m/s')

# === Set arrow legend value (uniformly set to 120) === 
arrow_value = 15  # Optional 100, 120, 130, uniform for dry and wet years for comparison

# === Plotting ===
ax_d = fig.add_subplot(2, 2, 4, projection=ccrs.PlateCarree())

# Read and draw shp
china = shpreader.Reader('/home/tangzhen/northwest/northwest_China.shp').geometries()
ax_d.add_geometries(china, ccrs.PlateCarree(), facecolor='none', edgecolor='green', linewidth=1.5, zorder=2, linestyle='-')

# Add geographic features
ax_d.coastlines()
ax_d.add_feature(cfeature.BORDERS, linestyle=':')
ax_d.add_feature(cfeature.LAND, facecolor='none')
ax_d.add_feature(cfeature.OCEAN, facecolor='white')
ax_d.set_extent([60, 120, 20, 60], crs=ccrs.PlateCarree())

# Set longitude and latitude ticks
xticks = np.arange(0, 181, 20)  # Remove 60 and 120
yticks = np.arange(0, 81, 20)   # Remove 20 and 60
ax_d.set_xticks(xticks, crs=ccrs.PlateCarree())
ax_d.set_yticks(yticks, crs=ccrs.PlateCarree())
ax_d.set_xticks(xticks, crs=ccrs.PlateCarree())
ax_d.set_yticks(yticks, crs=ccrs.PlateCarree())

# Set tick label format
ax_d.set_xticklabels([f'{lon}°E' for lon in xticks], fontsize=18)
ax_d.set_yticklabels([f'{lat}°N' for lat in yticks], fontsize=18)

# Show tick lines (ticks)
ax_d.tick_params(axis='both', which='major', length=5, width=1, labelsize=10)

# === Draw uniformly scaled vector arrows ===
q = ax_d.quiver(
    lon_sparse, lat_sparse, u_sparse, v_sparse,
    transform=ccrs.PlateCarree(),
    width=0.002,
    regrid_shape=20,
    scale=400# Uniform scaling ratio
)

# === Add arrow reference legend ===
ax_d.quiverkey(
    q, 0.9, 1.02, arrow_value,
    f'{arrow_value}m/s',
    labelpos='E', coordinates='axes',
    fontproperties={'size': 12}
)

# === Add 12500 geopotential height contour ===
contour_levels = [5880]  # South Asian high corresponds to geopotential height 12500
cs = ax_d.contour(
    ds_z_subset.lon, ds_z_subset.lat, z_200hpa, levels=contour_levels,
    colors='red', linewidths=2, transform=ccrs.PlateCarree()
)
ax_d.clabel(cs, inline=True, fontsize=10, fmt='%d')

# === Add title and color bar ===
ax_d.set_title(title, fontsize=18, loc='center')
ax_d.text(0.01, 1.01, "(d)", transform=plt.gca().transAxes, fontsize=18, ha='left', va='bottom', fontweight='bold')
# === Save figure ===
plt.savefig('', dpi=800, bbox_inches="tight")