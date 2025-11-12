import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import cmaps
from matplotlib.colors import BoundaryNorm

lon_slice = slice(60, 120)
lat_slice = slice(30, 50) 
fig = plt.figure(figsize=(30, 10))

# Read topographic data
t_file = '/home/tangzhen/geo.nc'
ts =xr.open_dataset(t_file)
t = ts['z']
t = t.reindex(latitude=t.latitude[::-1])  # Reverse latitude
# Get longitude from t dataset
lon = t['longitude'].values  # Assume longitude variable name is 'longitude'
# Convert longitude from 0 to 259 to -180 to 180
lon = np.where(lon > 180, lon - 360, lon)
# Update longitude coordinates in t dataset
t['longitude'] = ('longitude', lon)
t_subset = t.sel(longitude=lon_slice)
t_subset = t_subset.sel(latitude=lat_slice)
pressures = t_subset / 9.8
pressures = 1013 * (1 - 6.5 / 288000 * pressures) ** 5.255
print(pressures)
pressures_date = pressures.max(dim='latitude')
pressures_date = pressures_date.squeeze()  # Remove all dimensions with length 1

# ------------------------ (a)----------------------
# Read data
file_path = '/data1/omage_mean/wet_omage_mean.nc'  # Replace with your file path
ds = xr.open_dataset(file_path)

# Extract required variables
w = ds['w']  # Vertical velocity variable
u = ds['u']  # Zonal wind speed variable
v = ds['v']
lon = ds['lon']
lat = ds['lat']
lev = ds['lev']  # Pressure levels

# Select longitude range from 60°E to 120°E
lon_slice = slice(60, 120)  # Set longitude range
w_subset = w.sel(lon=lon_slice)
u_subset = u.sel(lon=lon_slice)
v_subset = v.sel(lon=lon_slice)

# Reverse pressure level data so 1000 hPa is at bottom
lev_reversed = lev.values[::-1]

# Limit pressure level range from 1000 hPa to 100 hPa
pressure_range = lev_reversed[(lev_reversed >= 100) & (lev_reversed <= 1000)]

# Select zonal cross-section area (assume selecting a latitude range)
lat_slice = slice(30, 50)  # Adjust latitude range as needed
w_subset = w_subset.sel(lat=lat_slice)
u_subset = u_subset.sel(lat=lat_slice)
v_subset = v_subset.sel(lat=lat_slice)

# Reverse vertical velocity data to ensure consistency with reversed pressure level data
w_reversed = w_subset.sel(lev=pressure_range)
u_reversed = u_subset.sel(lev=pressure_range)
v_reversed = v_subset.sel(lev=pressure_range)

magnitude = np.sqrt(u_reversed**2 + v_reversed**2)
mean_flux = float(magnitude.mean())
print(f'Average water vapor flux intensity: {mean_flux:.2f} kg/m/s')
arrow_value = 8

# === Calculate maximum and minimum values of data ===
w_max = np.nanmax(np.abs(w_reversed.values))  # Use np.nanmax to exclude NaN
w_min = np.nanmin(np.abs(w_reversed.values))  # Use np.nanmin to exclude NaN


# === Output vmin and vmax, ensure they are scalars ===
print(f"vmin: {w_min}, vmax: {w_max}")

# === Output vmin and vmax, ensure they are scalars ===
print(f"vmin: {w_min}, vmax: {w_max}")


# Create figure
ax_a = fig.add_subplot(1, 2, 1)
# Set missing value areas to black
# Set missing value background color to white
ax_a.set_facecolor('white')

# === Set color bar breakpoints ===
boundaries1 = np.arange(-0.10, 0, 0.02)
boundaries2 = np.arange(0, 0.11 ,0.02)
boundaries = np.concatenate((boundaries1, boundaries2))
norm = BoundaryNorm(boundaries, ncolors=64)
# Draw zonal cross-section using reversed pressure level data
a = ax_a.contourf(w_subset.lon, pressure_range, w_reversed.mean(dim='lat',skipna=True), 
                cmap=cmaps.cmp_b2r,
                levels=boundaries,norm=norm,
                extend ='both'
                )
# Set labels and title
ax_a.set_xlabel('Longitude',fontsize=20)
ax_a.set_ylabel('Pressure Level (hPa)', fontsize=20)
ax_a.set_title('High-variability years',loc='center', fontsize=22)

# Set vertical axis range from 1000 hPa to 100 hPa
ax_a.set_ylim(1000, 100)  # Set range between 1000 hPa and 100 hPa

# Set vertical axis ticks to show required pressure levels
ax_a.set_yticks([1000, 850, 700, 500, 300, 200, 100])  # Can adjust pressure levels as needed

# Set longitude range from 60°E to 120°E
xticks = np.arange(60, 160, 10)  # Remove 60 and 120
# Set x and y axis tick label font sizes
ax_a.tick_params(axis='x', labelsize=16)  # Set x-axis tick label font size to 14
ax_a.tick_params(axis='y', labelsize=16)  # Set y-axis tick label font size to 14
ax_a.set_xticklabels([f'{lon}°E' for lon in xticks], fontsize=16)

# Create grid
lon_grid, lev_grid = np.meshgrid(u_subset.lon.values, pressure_range)
# Add wind field arrows - only sparse longitude
quiver_skip = 8
q_a = ax_a.quiver(
    lon_grid[:, ::quiver_skip], 
    lev_grid[:, ::quiver_skip],
    u_reversed.mean(dim='lat').values[:, ::quiver_skip],
    v_reversed.mean(dim='lat').values[:, ::quiver_skip],
    scale=170, 
    color='black', 
    width=0.003
)
# === Add arrow reference legend ===
ax_a.quiverkey(
    q_a, 1.1, 0.02, arrow_value,
    f'{arrow_value} m/s',
    labelpos='E', coordinates='axes',
    fontproperties={'size': 16}
)

ax_a.fill_between(pressures.longitude, pressures_date ,1000,where=pressures_date<1010,facecolor='k') 
# Add color bar
# Add color bar
cb_a = fig.colorbar(a, ax=ax_a, label='Vertical Velocity (m/s)', orientation='vertical',shrink=0.8,aspect=20)

# Adjust color bar tick font size
cb_a.ax.tick_params(labelsize=16)

# If you want to adjust color bar label font size, you can do this:
cb_a.set_label('Vertical Velocity (m/s)', fontsize=20)
ax_a.text(0.00, 1.01, "(a)", transform=ax_a.transAxes, fontsize=22, ha='left', va='bottom', fontweight='bold')

# ------------------------ (b)----------------------
# Read data
file_path = '/data1/omage_mean/dry_omage_mean.nc'  # Replace with your file path
ds = xr.open_dataset(file_path)

# Extract required variables
w = ds['w']  # Vertical velocity variable
u = ds['u']  # Zonal wind speed variable
v = ds['v']
lon = ds['lon']
lat = ds['lat']
lev = ds['lev']  # Pressure levels

# Select longitude range from 60°E to 120°E
lon_slice = slice(60, 120)  # Set longitude range
w_subset = w.sel(lon=lon_slice)
u_subset = u.sel(lon=lon_slice)
v_subset = v.sel(lon=lon_slice)

# Reverse pressure level data so 1000 hPa is at bottom
lev_reversed = lev.values[::-1]

# Limit pressure level range from 1000 hPa to 100 hPa
pressure_range = lev_reversed[(lev_reversed >= 100) & (lev_reversed <= 1000)]

# Select zonal cross-section area (assume selecting a latitude range)
lat_slice = slice(30, 50)  # Adjust latitude range as needed
w_subset = w_subset.sel(lat=lat_slice)
u_subset = u_subset.sel(lat=lat_slice)
v_subset = v_subset.sel(lat=lat_slice)

# Reverse vertical velocity data to ensure consistency with reversed pressure level data
w_reversed = w_subset.sel(lev=pressure_range)
u_reversed = u_subset.sel(lev=pressure_range)
v_reversed = v_subset.sel(lev=pressure_range)

magnitude = np.sqrt(u_reversed**2 + v_reversed**2)
mean_flux = float(magnitude.mean())
print(f'Average water vapor flux intensity: {mean_flux:.2f} kg/m/s')
arrow_value = 8

# === Calculate maximum and minimum values of data ===
w_max = np.nanmax(np.abs(w_reversed.values))  # Use np.nanmax to exclude NaN
w_min = np.nanmin(np.abs(w_reversed.values))  # Use np.nanmin to exclude NaN


# === Output vmin and vmax, ensure they are scalars ===
print(f"vmin: {w_min}, vmax: {w_max}")

# === Output vmin and vmax, ensure they are scalars ===
print(f"vmin: {w_min}, vmax: {w_max}")


# Create figure
ax = fig.add_subplot(1, 2, 2)
ax.set_facecolor('white')

boundaries1 = np.arange(-0.10, 0, 0.02)
boundaries2 = np.arange(0, 0.11 ,0.02)
boundaries = np.concatenate((boundaries1, boundaries2))
norm = BoundaryNorm(boundaries, ncolors=64)
c = ax.contourf(w_subset.lon, pressure_range, w_reversed.mean(dim='lat', skipna=True), 
                cmap=cmaps.cmp_b2r,
                levels=boundaries,norm=norm,
                extend ='both'
                )
# Set labels and title
ax.set_xlabel('Longitude',fontsize=20)
ax.set_ylabel('Pressure Level (hPa)', fontsize=20)
ax.set_title('Low-variability years',loc='center', fontsize=22)

# Set vertical axis range from 1000 hPa to 100 hPa
ax.set_ylim(1000, 100)  # Set range between 1000 hPa and 100 hPa

# Set vertical axis ticks to show required pressure levels
ax.set_yticks([1000, 850, 700, 500, 300, 200, 100])  # Can adjust pressure levels as needed

# Set longitude range from 60°E to 120°E
xticks = np.arange(60, 160, 10)  # Remove 60 and 120
# Set x and y axis tick label font sizes
ax.tick_params(axis='x', labelsize=16)  # Set x-axis tick label font size to 14
ax.tick_params(axis='y', labelsize=16)  # Set y-axis tick label font size to 14
ax.set_xticklabels([f'{lon}°E' for lon in xticks], fontsize=16)

# Create grid
lon_grid, lev_grid = np.meshgrid(u_subset.lon.values, pressure_range)
# Add wind field arrows - only sparse longitude
quiver_skip = 8
q = ax.quiver(
    lon_grid[:, ::quiver_skip], 
    lev_grid[:, ::quiver_skip],
    u_reversed.mean(dim='lat').values[:, ::quiver_skip],
    v_reversed.mean(dim='lat').values[:, ::quiver_skip],
    scale=170, 
    color='black', 
    width=0.003
)
# === Add arrow reference legend ===
ax.quiverkey(
    q, 1.1, 0.02, arrow_value,
    f'{arrow_value} m/s',
    labelpos='E', coordinates='axes',
    fontproperties={'size': 16}
)

ax.fill_between(pressures.longitude, pressures_date ,1000,where=pressures_date<1010,facecolor='k') 

cbar = fig.colorbar(c, ax=ax, label='Vertical Velocity (m/s)', orientation='vertical',shrink=0.8,aspect=20)

# Adjust color bar tick font size
cbar.ax.tick_params(labelsize=16)

# If you want to adjust color bar label font size, you can do this:
cbar.set_label('Vertical Velocity (m/s)', fontsize=20)
ax.text(0.00, 1.01, "(b)", transform=plt.gca().transAxes, fontsize=22, ha='left', va='bottom', fontweight='bold')

plt.savefig('', dpi=600, bbox_inches="tight")