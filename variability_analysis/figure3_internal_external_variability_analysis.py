import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy.stats import linregress
import cartopy.io.shapereader as shpreader
from netCDF4 import Dataset
from matplotlib.colors import BoundaryNorm
import cartopy.crs as ccrs
import cartopy.mpl.ticker as cticker
import cmaps

# Set region boundaries
leftlon, rightlon, lowerlat, upperlat = 70, 115, 30, 51

# === Set up figure ===
fig = plt.figure(figsize=(24, 24))

#=== (a) internal_prev_plot ===
folder_path = ''
# Get paths for all .nc files
file_paths = [os.path.join(folder_path, f'internal_prev_{year}.nc') for year in range(1982, 2021)]
# Calculate weighted average for each year
pre_int_mean = []
for file_path in file_paths:
    ds = xr.open_dataset(file_path)
    data = ds['tp']*1000
    weights = np.cos(np.deg2rad(data.latitude))
    weitght_avg = data.weighted(weights).mean(dim=["latitude","longitude"])
    pre_int_mean.append(weitght_avg)

mean_value = np.mean(pre_int_mean)
print(mean_value)
# Calculate linear trend
years = np.arange(1982,2021)
slop, intercept, r_value, p_value, std_err = linregress(years, pre_int_mean)
trend_line = slop*years + intercept

ax_a = fig.add_subplot(3, 2, 1)

# Main plot content
ax_a.plot(years, pre_int_mean, marker='s', color='red', markersize=10, linewidth=2, linestyle='--')
ax_a.plot(years, trend_line, linestyle='-', color='black', linewidth=2)

# Axis settings
ax_a.xaxis.set_minor_locator(MultipleLocator(1))
ax_a.xaxis.set_major_locator(MultipleLocator(5))
ax_a.tick_params(axis='both', which='minor', width=1, length=4, labelsize=16)
ax_a.tick_params(axis='both', which='major', width=1, length=6, labelsize=16)

# Annotate trend
ax_a.annotate(f'trend = {slop:.4f} mm/year (p < 0.05)',
              xy=(years[20], trend_line[20]),
              xytext=(-75, -120),
              textcoords='offset points',
              fontsize=20,
              color='red'
              )

# Add title, axis labels, and panel label
ax_a.set_title('Internal', loc='center', fontsize=20)
ax_a.set_ylabel('mm', fontsize=18)
ax_a.set_xlabel('Years', fontsize=18)
ax_a.xaxis.set_label_coords(0.5, -0.12)  # Horizontal 0.5 center, vertical -0.1下移
ax_a.text(0.00, 1.02, "(a)", transform=ax_a.transAxes, fontsize=20, ha='left', va='bottom', fontweight='bold')

#=== (b) external prev plot ===
folder_path = ''
file_paths = [os.path.join(folder_path, f'prev_{year}.nc') for year in range(1982, 2021)]
pre_int_mean = []
for file_path in file_paths:
    ds = xr.open_dataset(file_path)
    data = ds['tp']
    weights = np.cos(np.deg2rad(data.latitude))
    weitght_avg = data.weighted(weights).mean(dim=["latitude","longitude"])
    pre_int_mean.append(weitght_avg)

years = np.arange(1982,2021)
slop, intercept, r_value, p_value, std_err = linregress(years, pre_int_mean)
trend_line = slop*years + intercept

ax_b = fig.add_subplot(3, 2, 2)
ax_b.plot(years, pre_int_mean, marker='s', color='blue', markersize=10, linewidth=2, linestyle='--')

ax_b.plot(years, trend_line, linestyle='-', color='black', linewidth=2)
ax_b.xaxis.set_minor_locator(MultipleLocator(1))
ax_b.xaxis.set_major_locator(MultipleLocator(5))
ax_a.tick_params(axis='both', which='minor', width=1, length=4, labelsize=16)
ax_a.tick_params(axis='both', which='major', width=1, length=6, labelsize=16)

ax_b.tick_params(axis='both',
                which='minor',
                width=1,
                length=4,
                labelsize=16)
ax_b.tick_params(axis='both',
                which='major',
                width=1,
                length=6,
                labelsize=18)
ax_b.annotate(f'trend={slop:.4f}mm/year(p<0.05)',
             xy=(years[20],trend_line[20]),
             xytext=(-50,-120),
             textcoords='offset points',
             color='blue',
             fontsize=20)
ax_b.set_title('External',loc='center',fontsize=20)
ax_b.set_ylabel('mm',fontsize=18)
ax_b.set_xlabel('Years',fontsize=18)
ax_b.xaxis.set_label_coords(0.5, -0.12)
ax_b.text(0.00, 1.02, "(b)", transform=plt.gca().transAxes, fontsize=20, ha='left', va='bottom', fontweight='bold')

#=== (c) internal prev contourf ===
file_path = ''
nc_data = Dataset(file_path, 'r')

var_name = 'tp'
variable = nc_data.variables[var_name][:]*1000

lon = nc_data.variables['longitude'][:]
lat = nc_data.variables['latitude'][:] 
ax_c = fig.add_subplot(3, 2, 3, projection=ccrs.PlateCarree(central_longitude=90))
ax_c.set_extent([leftlon, rightlon, lowerlat, upperlat], crs=ccrs.PlateCarree())

xticks = np.arange(leftlon, rightlon, 10)
yticks = np.arange(lowerlat, upperlat, 5)
xticks = xticks[1:]  # Exclude first longitude
yticks = yticks[1:]  # Exclude first latitude

ax_c.set_xticks(xticks, crs=ccrs.PlateCarree())
ax_c.set_yticks(yticks, crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax_c.xaxis.set_major_formatter(lon_formatter)
ax_c.yaxis.set_major_formatter(lat_formatter)
ax_c.tick_params(axis='x', labelsize=16)  # x-axis font size
ax_c.tick_params(axis='y', labelsize=16)  # y-axis font size
ax_c.set_title('Climatology of internal precipitation variability',loc='center',fontsize=20)
ax_c.text(0.00, 1.02, "(c)", transform=ax_c.transAxes, fontsize=20, ha='left', va='bottom', fontweight='bold')

china = shpreader.Reader('/home/tangzhen/northwest/northwest_China.shp').geometries()
ax_c.add_geometries(china, ccrs.PlateCarree(),facecolor='none', edgecolor='black', linewidth=1.5, zorder = 2)

boundaries = np.arange(0.0, 1.21, 0.1)
norm = BoundaryNorm(boundaries, ncolors=128)
# Plot spatial distribution of tp-intensity
internal_mean = ax_c.contourf(lon, lat,variable, cmap=cmaps.MPL_YlGnBu,extend='max', levels=boundaries,norm=norm, transform=ccrs.PlateCarree(),zorder=1)

cb_ax_c = fig.add_axes([0.13, 0.36, 0.33, 0.01])
cb_c = fig.colorbar(internal_mean, cax=cb_ax_c, orientation='horizontal')
cb_c.set_label('mm', fontsize=18)
cb_c.set_ticks(np.arange(0.0 ,1.21, 0.1))  # Set ticks every 1.0
cb_c.ax.tick_params(labelsize=14)  # Set colorbar tick label font size to 14 

# ===(d) external contourf ===
file_path = ''
nc_data = Dataset(file_path, 'r')

var_name = 'tp'
variable = nc_data.variables[var_name][:]

lon = nc_data.variables['longitude'][:]
lat = nc_data.variables['latitude'][:] 

ax_d = fig.add_subplot(3, 2, 4, projection=ccrs.PlateCarree(central_longitude=90))
ax_d.set_extent([leftlon, rightlon, lowerlat, upperlat], crs=ccrs.PlateCarree())
xticks = np.arange(leftlon, rightlon, 10)
yticks = np.arange(lowerlat, upperlat, 5)

xticks = xticks[1:]  
yticks = yticks[1:] 

ax_d.set_xticks(xticks, crs=ccrs.PlateCarree())
ax_d.set_yticks(yticks, crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax_d.xaxis.set_major_formatter(lon_formatter)
ax_d.yaxis.set_major_formatter(lat_formatter)
ax_d.tick_params(axis='x', labelsize=16)  # x-axis font size
ax_d.tick_params(axis='y', labelsize=16)  # y-axis font size
ax_d.set_title('Climatology of external precipitation variability',loc='center',fontsize=20)
ax_d.text(0.00, 1.02, "(d)", transform=plt.gca().transAxes, fontsize=20, ha='left', va='bottom', fontweight='bold')

china = shpreader.Reader('/home/tangzhen/northwest/northwest_China.shp').geometries()
ax_d.add_geometries(china, ccrs.PlateCarree(),facecolor='none', edgecolor='black', linewidth=1.5, zorder = 2)

boundaries = np.arange(0.0, 8.1, 0.5)
norm = BoundaryNorm(boundaries, ncolors=128)
tp_int = ax_d.contourf(lon, lat,variable, cmap=cmaps.MPL_YlGnBu,extend='max', levels=boundaries,norm=norm, transform=ccrs.PlateCarree(),zorder=1)
cb_ax_d = fig.add_axes([0.57, 0.36, 0.33, 0.01])
cb_d = fig.colorbar(tp_int, cax=cb_ax_d, orientation='horizontal')
cb_d.set_label('mm', fontsize=18)
cb_d.set_ticks(np.arange(0.0, 8.1, 0.5))  # Set ticks every 1.0
cb_d.ax.tick_params(labelsize=14)  # Set colorbar tick label font size to 14

# ===(e) internal_trend
file_path = ''
nc_data = Dataset(file_path, 'r')

var_name = 'trend'
variable = nc_data.variables[var_name][:]*100000

p_value = nc_data.variables['p_value'][:]
lon = nc_data.variables['longitude'][:]
lat = nc_data.variables['latitude'][:] 

lon_grid, lat_grid = np.meshgrid(lon, lat)

ax_e = fig.add_subplot(3, 2, 5, projection=ccrs.PlateCarree(central_longitude=90))
ax_e.set_extent([leftlon, rightlon, lowerlat, upperlat], crs=ccrs.PlateCarree())
xticks = np.arange(leftlon, rightlon, 10)
yticks = np.arange(lowerlat, upperlat, 5)

xticks = xticks[1:]
yticks = yticks[1:] 

ax_e.set_xticks(xticks, crs=ccrs.PlateCarree())
ax_e.set_yticks(yticks, crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax_e.xaxis.set_major_formatter(lon_formatter)
ax_e.yaxis.set_major_formatter(lat_formatter)
ax_e.tick_params(axis='x', labelsize=16)  # x-axis font size
ax_e.tick_params(axis='y', labelsize=16)  # y-axis font size
ax_e.set_title('Trend of internal precipitation variability',loc='center',fontsize=20)
ax_e.text(0.00, 1.02, "(e)", transform=plt.gca().transAxes, fontsize=20, ha='left', va='bottom', fontweight='bold')

china = shpreader.Reader('/home/northwest/northwest_China.shp').geometries()
ax_e.add_geometries(china, ccrs.PlateCarree(),facecolor='none', edgecolor='black', linewidth=1.5, zorder = 2)


boundaries1 = np.arange(-0.05, 0, 0.005)
boundaries2 = np.arange(0, 1.5 ,0.15)
boundaries = np.concatenate((boundaries1, boundaries2))
norm = BoundaryNorm(boundaries, ncolors=128)
tp_int = ax_e.contourf(lon, lat,variable, cmap=cmaps.MPL_BrBG,extend='both', levels=boundaries,norm=norm, transform=ccrs.PlateCarree(),zorder=1)

significant_mask = p_value < 0.05
significant_lons = lon_grid[significant_mask]
significant_lats = lat_grid[significant_mask]
ax_e.scatter(significant_lons, significant_lats, color='black', s=0.05, transform=ccrs.PlateCarree(), zorder=3,marker='x')


cb_ax_e = fig.add_axes([0.13, 0.09, 0.33, 0.01])
cb_e = fig.colorbar(tp_int, cax=cb_ax_e, orientation='horizontal')
cb_e.set_label('mm/year (%)', fontsize=18)
cb_e.ax.tick_params(labelsize=14)  # Set colorbar tick label font size to 14

#=== external trend ===
# Read nc file
file_path = ''
nc_data = Dataset(file_path, 'r')

# Plot spatial distribution of variable tp
var_name = 'trend'
variable = nc_data.variables[var_name][:]*100

print("Min:", np.min(variable))
print("Max:", np.max(variable))


# Get p_value data
p_value = nc_data.variables['p_value'][:]

# Get longitude and latitude
lon = nc_data.variables['longitude'][:]
lat = nc_data.variables['latitude'][:] 

# Create longitude-latitude grid
lon_grid, lat_grid = np.meshgrid(lon, lat)

# Draw map
ax_f = fig.add_subplot(3, 2, 6, projection=ccrs.PlateCarree(central_longitude=90))
ax_f.set_extent([leftlon, rightlon, lowerlat, upperlat], crs=ccrs.PlateCarree())
xticks = np.arange(leftlon, rightlon, 10)
yticks = np.arange(lowerlat, upperlat, 5)

# Remove first coordinate point
xticks = xticks[1:]  # Exclude first longitude
yticks = yticks[1:]  # Exclude first latitude

ax_f.set_xticks(xticks, crs=ccrs.PlateCarree())
ax_f.set_yticks(yticks, crs=ccrs.PlateCarree())
# Format longitude and latitude tick labels to standard format
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax_f.xaxis.set_major_formatter(lon_formatter)
ax_f.yaxis.set_major_formatter(lat_formatter)
ax_f.tick_params(axis='x', labelsize=16)  # x-axis font size
ax_f.tick_params(axis='y', labelsize=16)  # y-axis font size
ax_f.set_title('Trend of external precipitation variability',loc='center',fontsize=18)
ax_f.text(0.00, 1.02, "(f)", transform=plt.gca().transAxes, fontsize=20, ha='left', va='bottom', fontweight='bold')

china = shpreader.Reader('/home/tangzhen/northwest/northwest_China.shp').geometries()
# Draw China border and nine-dash line
ax_f.add_geometries(china, ccrs.PlateCarree(),facecolor='none', edgecolor='black', linewidth=1.5, zorder = 2)

# Set colorbar breakpoints
boundaries1 = np.arange(-1.5, 0, 0.15)
boundaries2 = np.arange(0, 5.1 ,0.5)
boundaries = np.concatenate((boundaries1, boundaries2))
norm = BoundaryNorm(boundaries, ncolors=128)
# Plot spatial distribution of tp-intensity
tp_int = ax_f.contourf(lon, lat,variable, cmap=cmaps.MPL_BrBG,extend='both', levels=boundaries,norm=norm, transform=ccrs.PlateCarree(),zorder=1)
#tp_int = f2_ax1.contourf(lon, lat,variable, cmap=cmaps.MPL_BrBG,extend='both',  transform=ccrs.PlateCarree())

# Mark grid points with p_value less than 0.05
significant_mask = p_value < 0.05
# Get longitude and latitude of significant grid points
significant_lons = lon_grid[significant_mask]
significant_lats = lat_grid[significant_mask]
# Mark significant grid points with small black dots on the map
ax_f.scatter(significant_lons, significant_lats, color='black', s=0.05, transform=ccrs.PlateCarree(), zorder=3,marker='x')


# Add colorbar
cb_ax_f = fig.add_axes([0.57, 0.09, 0.33, 0.01])
cb_f = fig.colorbar(tp_int, cax=cb_ax_f, orientation='horizontal')
cb_f.set_label('mm/year (%)', fontsize=18)
cb_f.ax.tick_params(labelsize=14)  # Set colorbar tick label font size to 14



plt.subplots_adjust(hspace=0.25, wspace=0.25)
plt.savefig('', dpi=800, bbox_inches='tight')
plt.show()