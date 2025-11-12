import numpy as np
import cmaps
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.mpl.ticker as cticker
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature
from netCDF4 import Dataset
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

leftlon, rightlon, lowerlat, upperlat = 70, 115, 30, 51
fig = plt.figure(figsize=(24, 30))

# ------------------------ (a) NW ------------------------
file_path = '/data/nw_prev/mean_prev_1982_2020.nc'
nc_data = Dataset(file_path, 'r')
var_name = 'tp'
variable = nc_data.variables[var_name][:]
print("Min:", np.min(variable))
print("Max:", np.max(variable))


lon = nc_data.variables['longitude'][:]
lat = nc_data.variables['latitude'][:] 
ax_a = fig.add_subplot(4, 2, 1, projection=ccrs.PlateCarree())
ax_a.set_extent([leftlon, rightlon, lowerlat, upperlat], crs=ccrs.PlateCarree())
xticks = np.arange(leftlon, rightlon, 10)
yticks = np.arange(lowerlat, upperlat, 5)
xticks = xticks[1:]  
yticks = yticks[1:] 

ax_a.set_xticks(xticks, crs=ccrs.PlateCarree())
ax_a.set_yticks(yticks, crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax_a.xaxis.set_major_formatter(lon_formatter)
ax_a.yaxis.set_major_formatter(lat_formatter)
ax_a.tick_params(axis='x', labelsize=16)  # X-axis font size
ax_a.tick_params(axis='y', labelsize=16)  # Y-axis font size
ax_a.set_title('Climatological of precipitation variability in the NW', loc='center', fontsize=20)

china = shpreader.Reader('/home/northwest/northwest_China.shp').geometries()
ax_a.add_geometries(china, ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=1.5, zorder=2)

boundaries = np.arange(0.2, 2.1, 0.2)
norm = BoundaryNorm(boundaries, ncolors=128)
tp_int = ax_a.contourf(lon, lat, variable, cmap=cmaps.MPL_YlGnBu, extend='both', levels=boundaries, norm=norm, transform=ccrs.PlateCarree(), zorder=1)
cb_ax_a = fig.add_axes([0.13, 0.70, 0.33, 0.01])
cb_a = fig.colorbar(tp_int, cax=cb_ax_a, orientation='horizontal')
cb_a.set_label('mm', fontsize=16)
cb_a.ax.tick_params(labelsize=16)  # Set colorbar tick label font size to 12
#ax_a.text(-0.05, 16.30, "(a)", transform=plt.gca().transAxes, fontsize=18, ha='left', va='bottom', fontweight='bold')

# ------------------------ (b) NE ------------------------
file_path = '/data/ne_prev/mean_prev_1982_2020.nc'
nc_data = Dataset(file_path, 'r')
var_name = 'tp'
variable = nc_data.variables[var_name][:]
print("Min:", np.min(variable))
print("Max:", np.max(variable))

lon = nc_data.variables['longitude'][:]
lat = nc_data.variables['latitude'][:] 

ax_b = fig.add_subplot(4, 2, 2, projection=ccrs.PlateCarree())
ax_b.set_extent([leftlon, rightlon, lowerlat, upperlat], crs=ccrs.PlateCarree())
xticks = np.arange(leftlon, rightlon, 10)
yticks = np.arange(lowerlat, upperlat, 5)

xticks = xticks[1:] 
yticks = yticks[1:]  

ax_b.set_xticks(xticks, crs=ccrs.PlateCarree())
ax_b.set_yticks(yticks, crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax_b.xaxis.set_major_formatter(lon_formatter)
ax_b.yaxis.set_major_formatter(lat_formatter)
ax_b.tick_params(axis='x', labelsize=16)  # X-axis font size
ax_b.tick_params(axis='y', labelsize=16)  # Y-axis font size
ax_b.set_title('Climatological of precipitation variability in the NE', loc='center', fontsize=20)

china = shpreader.Reader('/home/tangzhen/northwest/northwest_China.shp').geometries()
ax_b.add_geometries(china, ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=1.5, zorder=2)

boundaries = np.arange(0.05, 0.51, 0.05)
norm = BoundaryNorm(boundaries, ncolors=128)
tp_int = ax_b.contourf(lon, lat, variable, cmap=cmaps.MPL_YlGnBu, extend='both', levels=boundaries, norm=norm, transform=ccrs.PlateCarree(), zorder=1)
cb_ax_b = fig.add_axes([0.56, 0.70, 0.33, 0.01])
cb_b = fig.colorbar(tp_int, cax=cb_ax_b, orientation='horizontal')
cb_b.set_label('mm', fontsize=16)
cb_b.ax.tick_params(labelsize=16)  # Set colorbar tick label font size to 12
#ax_b.text(0.00, 16.30, "(b)", transform=plt.gca().transAxes, fontsize=18, ha='left', va='bottom', fontweight='bold')


# ------------------------ (c) SW ------------------------
file_path = '/data/sw_prev/mean_prev_1982_2020.nc'
nc_data = Dataset(file_path, 'r')
var_name = 'tp'
variable = nc_data.variables[var_name][:]
print("Min:", np.min(variable))
print("Max:", np.max(variable))

lon = nc_data.variables['longitude'][:]
lat = nc_data.variables['latitude'][:] 
ax_c = fig.add_subplot(4, 2, 3, projection=ccrs.PlateCarree())
ax_c.set_extent([leftlon, rightlon, lowerlat, upperlat], crs=ccrs.PlateCarree())

#Define geographic coordinate label format
xticks = np.arange(leftlon, rightlon, 10)
yticks = np.arange(lowerlat, upperlat, 5)

# Remove first coordinate point
xticks = xticks[1:]  # Exclude first longitude
yticks = yticks[1:]  # Exclude first latitude

ax_c.set_xticks(xticks, crs=ccrs.PlateCarree())
ax_c.set_yticks(yticks, crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax_c.xaxis.set_major_formatter(lon_formatter)
ax_c.yaxis.set_major_formatter(lat_formatter)
ax_c.tick_params(axis='x', labelsize=16)  # X-axis font size
ax_c.tick_params(axis='y', labelsize=16)  # Y-axis font size
ax_c.set_title('Climatological of precipitation variability in the SW', loc='center', fontsize=20)

china = shpreader.Reader('/home/tangzhen/northwest/northwest_China.shp').geometries()
ax_c.add_geometries(china, ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=1.5, zorder=2)

boundaries = np.arange(0.25, 5.1, 0.25)
norm = BoundaryNorm(boundaries, ncolors=128)
tp_int = ax_c.contourf(lon, lat, variable, cmap=cmaps.MPL_YlGnBu, extend='both', levels=boundaries, norm=norm, transform=ccrs.PlateCarree(), zorder=1)
cb_ax_c = fig.add_axes([0.13, 0.50, 0.33, 0.01])
cb_c = fig.colorbar(tp_int, cax=cb_ax_c, orientation='horizontal')
cb_c.set_label('mm', fontsize=16)
cb_c.ax.tick_params(labelsize=16)  # Set colorbar tick label font size to 12
#ax_c.text(-0.05, 16.30, "(c)", transform=plt.gca().transAxes, fontsize=18, ha='left', va='bottom', fontweight='bold')

# ------------------------ (d) MERRA2 Spatial Trend ------------------------
file_path = '/data/_se_prev/mean_prev_1982_2020.nc'
nc_data = Dataset(file_path, 'r')
var_name = 'tp'
variable = nc_data.variables[var_name][:]
print("Min:", np.min(variable))
print("Max:", np.max(variable))

lon = nc_data.variables['longitude'][:]
lat = nc_data.variables['latitude'][:] 
ax_d = fig.add_subplot(4, 2, 4, projection=ccrs.PlateCarree())
ax_d.set_extent([leftlon, rightlon, lowerlat, upperlat], crs=ccrs.PlateCarree())
xticks = np.arange(leftlon, rightlon, 10)
yticks = np.arange(lowerlat, upperlat, 5)
xticks = xticks[1:]  # Exclude first longitude
yticks = yticks[1:]  # Exclude first latitude

ax_d.set_xticks(xticks, crs=ccrs.PlateCarree())
ax_d.set_yticks(yticks, crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax_d.xaxis.set_major_formatter(lon_formatter)
ax_d.yaxis.set_major_formatter(lat_formatter)
ax_d.tick_params(axis='x', labelsize=16)  # X-axis font size
ax_d.tick_params(axis='y', labelsize=16)  # Y-axis font size
ax_d.set_title('Climatological of precipitation variability in the SE', loc='center', fontsize=20)

china = shpreader.Reader('/home/tangzhen/northwest/northwest_China.shp').geometries()
ax_d.add_geometries(china, ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=1.5, zorder=2)

boundaries = np.arange(0.1, 1.51, 0.1)
norm = BoundaryNorm(boundaries, ncolors=128)
tp_int = ax_d.contourf(lon, lat, variable, cmap=cmaps.MPL_YlGnBu, extend='both', levels=boundaries, norm=norm, transform=ccrs.PlateCarree(), zorder=1)
cb_ax_d = fig.add_axes([0.56, 0.50, 0.33, 0.01])
cb_d = fig.colorbar(tp_int, cax=cb_ax_d, orientation='horizontal')
cb_d.set_label('mm', fontsize=16)
cb_d.ax.tick_params(labelsize=16)  # Set colorbar tick label font size to 12
#ax_d.text(0.00, 16.30, "(d)", transform=plt.gca().transAxes, fontsize=18, ha='left', va='bottom', fontweight='bold')

# ------------------------ (e) MERRA2 Spatial Trend ------------------------
# Read nc file
file_path = '/data/nw_prev/prev_trend_and_p_1982_2020.nc'
nc_data = Dataset(file_path, 'r')

# Draw spatial distribution map of variable tp
var_name = 'trend'
variable = nc_data.variables[var_name][:]*100

print("Min:", np.min(variable))
print("Max:", np.max(variable))


# Get p_value data
p_value = nc_data.variables['p_value'][:]

# Get longitude and latitude
lon = nc_data.variables['longitude'][:]
lat = nc_data.variables['latitude'][:] 

# Create longitude and latitude grid
lon_grid, lat_grid = np.meshgrid(lon, lat)

# Draw map
ax_e = fig.add_subplot(4, 2, 5, projection=ccrs.PlateCarree())
# Create subplot in absolute coordinates of canvas
ax_e.set_extent([leftlon, rightlon, lowerlat, upperlat], crs=ccrs.PlateCarree())
# Define geographic coordinate label format
xticks = np.arange(leftlon, rightlon, 10)
yticks = np.arange(lowerlat, upperlat, 5)

# Remove first coordinate point
xticks = xticks[1:]  # Exclude first longitude
yticks = yticks[1:]  # Exclude first latitude

ax_e.set_xticks(xticks, crs=ccrs.PlateCarree())
ax_e.set_yticks(yticks, crs=ccrs.PlateCarree())
# Format longitude and latitude tick labels to standard format
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax_e.xaxis.set_major_formatter(lon_formatter)
ax_e.yaxis.set_major_formatter(lat_formatter)
ax_e.tick_params(axis='x', labelsize=16)  # X-axis font size
ax_e.tick_params(axis='y', labelsize=16)  # Y-axis font size
ax_e.set_title('Trend of precipitation variability in the NW',loc='center',fontsize=18)

# Read shp file
china = shpreader.Reader('/home/tangzhen/northwest/northwest_China.shp').geometries()
# Draw China national boundary with nine-dash line
ax_e.add_geometries(china, ccrs.PlateCarree(),facecolor='none', edgecolor='black', linewidth=1.5, zorder = 2)

# Set color bar breakpoints
boundaries1 = np.arange(-1.5, 0, 0.25)
boundaries2 = np.arange(0, 2.1 ,0.4)
boundaries = np.concatenate((boundaries1, boundaries2))
norm = BoundaryNorm(boundaries, ncolors=128)
# Draw spatial distribution map of tp-intensity
tp_int = ax_e.contourf(lon, lat,variable, cmap=cmaps.MPL_BrBG,extend='both', levels=boundaries,norm=norm, transform=ccrs.PlateCarree(),zorder=1)
#tp_int = f2_ax1.contourf(lon, lat,variable, cmap=cmaps.MPL_BrBG,extend='both',  transform=ccrs.PlateCarree())

# Mark grid points with p_value less than 0.05
significant_mask = p_value < 0.05
# Get longitude and latitude of significant grid points
significant_lons = lon_grid[significant_mask]
significant_lats = lat_grid[significant_mask]
# Mark significant grid points with small black dots on the map
ax_e.scatter(significant_lons, significant_lats, color='black', s=0.05, transform=ccrs.PlateCarree(), zorder=3,marker='x')
# Add color bar
cb_ax_e = fig.add_axes([0.13, 0.30, 0.33, 0.01])
cb_e = fig.colorbar(tp_int, cax=cb_ax_e, orientation='horizontal')
cb_e.set_label('mm/year (%)', fontsize=16)
cb_e.ax.tick_params(labelsize=12)  # Set colorbar tick label font size to 14
# ------------------------ (f) MERRA2 Spatial Trend ------------------------
# Read nc file
file_path = '/data/ne_prev/prev_trend_and_p_1982_2020.nc'
nc_data = Dataset(file_path, 'r')

# Draw spatial distribution map of variable tp
var_name = 'trend'
variable = nc_data.variables[var_name][:]*100

print("Min:", np.min(variable))
print("Max:", np.max(variable))


# Get p_value data
p_value = nc_data.variables['p_value'][:]

# Get longitude and latitude
lon = nc_data.variables['longitude'][:]
lat = nc_data.variables['latitude'][:] 

# Create longitude and latitude grid
lon_grid, lat_grid = np.meshgrid(lon, lat)

# Draw map
ax_f = fig.add_subplot(4, 2, 6, projection=ccrs.PlateCarree())
ax_f.set_extent([leftlon, rightlon, lowerlat, upperlat], crs=ccrs.PlateCarree())
# Define geographic coordinate label format
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
ax_f.tick_params(axis='x', labelsize=16)  # X-axis font size
ax_f.tick_params(axis='y', labelsize=16)  # Y-axis font size
ax_f.set_title('Trend of precipitation variability in the NE',loc='center',fontsize=18)

# Read shp file
china = shpreader.Reader('/home/tangzhen/northwest/northwest_China.shp').geometries()
# Draw China national boundary with nine-dash line
ax_f.add_geometries(china, ccrs.PlateCarree(),facecolor='none', edgecolor='black', linewidth=1.5, zorder = 2)

# Set color bar breakpoints
boundaries1 = np.arange(-0.4, 0, 0.04)
boundaries2 = np.arange(0, 1.1 ,0.1)
boundaries = np.concatenate((boundaries1, boundaries2))
norm = BoundaryNorm(boundaries, ncolors=128)
# Draw spatial distribution map of tp-intensity
tp_int = ax_f.contourf(lon, lat,variable, cmap=cmaps.MPL_BrBG,extend='both', levels=boundaries,norm=norm, transform=ccrs.PlateCarree(),zorder=1)
#tp_int = f2_ax1.contourf(lon, lat,variable, cmap=cmaps.MPL_BrBG,extend='both',  transform=ccrs.PlateCarree())

# Mark grid points with p_value less than 0.05
significant_mask = p_value < 0.05
# Get longitude and latitude of significant grid points
significant_lons = lon_grid[significant_mask]
significant_lats = lat_grid[significant_mask]
# Mark significant grid points with small black dots on the map
ax_f.scatter(significant_lons, significant_lats, color='black', s=0.05, transform=ccrs.PlateCarree(), zorder=3,marker='x')


# Add color bar
cb_ax_f = fig.add_axes([0.56, 0.30, 0.33, 0.01])
cb_f = fig.colorbar(tp_int, cax=cb_ax_f, orientation='horizontal')
cb_f.set_label('mm/year (%)', fontsize=16)
cb_f.ax.tick_params(labelsize=12)  # Set colorbar tick label font size to 14
# ------------------------ (g) MERRA2 Spatial Trend ------------------------
# Read nc file
file_path = '/data/sw_prev/prev_trend_and_p_1982_2020.nc'
nc_data = Dataset(file_path, 'r')

# Draw spatial distribution map of variable tp
var_name = 'trend'
variable = nc_data.variables[var_name][:]*100

print("Min:", np.min(variable))
print("Max:", np.max(variable))


# Get p_value data
p_value = nc_data.variables['p_value'][:]

# Get longitude and latitude
lon = nc_data.variables['longitude'][:]
lat = nc_data.variables['latitude'][:] 

# Create longitude and latitude grid
lon_grid, lat_grid = np.meshgrid(lon, lat)

# Draw map
ax_g = fig.add_subplot(4, 2, 7, projection=ccrs.PlateCarree())
ax_g.set_extent([leftlon, rightlon, lowerlat, upperlat], crs=ccrs.PlateCarree())
# Define geographic coordinate label format
xticks = np.arange(leftlon, rightlon, 10)
yticks = np.arange(lowerlat, upperlat, 5)

# Remove first coordinate point
xticks = xticks[1:]  # Exclude first longitude
yticks = yticks[1:]  # Exclude first latitude

ax_g.set_xticks(xticks, crs=ccrs.PlateCarree())
ax_g.set_yticks(yticks, crs=ccrs.PlateCarree())
# Format longitude and latitude tick labels to standard format
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax_g.xaxis.set_major_formatter(lon_formatter)
ax_g.yaxis.set_major_formatter(lat_formatter)
ax_g.tick_params(axis='x', labelsize=16)  # X-axis font size
ax_g.tick_params(axis='y', labelsize=16)  # Y-axis font size
ax_g.set_title('Trend of precipitation variability in the SW',loc='center',fontsize=18)

# Read shp file
china = shpreader.Reader('/home/tangzhen/northwest/northwest_China.shp').geometries()
# Draw China national boundary with nine-dash line
ax_g.add_geometries(china, ccrs.PlateCarree(),facecolor='none', edgecolor='black', linewidth=1.5, zorder = 2)

# Set color bar breakpoints
boundaries1 = np.arange(-1.5, 0, 0.15)
boundaries2 = np.arange(0, 3.61 ,0.36)
boundaries = np.concatenate((boundaries1, boundaries2))
norm = BoundaryNorm(boundaries, ncolors=128)
# Draw spatial distribution map of tp-intensity
tp_int = ax_g.contourf(lon, lat,variable, cmap=cmaps.MPL_BrBG,extend='both', levels=boundaries,norm=norm, transform=ccrs.PlateCarree(),zorder=1)
#tp_int = f2_ax1.contourf(lon, lat,variable, cmap=cmaps.MPL_BrBG,extend='both',  transform=ccrs.PlateCarree())

# Mark grid points with p_value less than 0.05
significant_mask = p_value < 0.05
# Get longitude and latitude of significant grid points
significant_lons = lon_grid[significant_mask]
significant_lats = lat_grid[significant_mask]
# Mark significant grid points with small black dots on the map
ax_g.scatter(significant_lons, significant_lats, color='black', s=0.05, transform=ccrs.PlateCarree(), zorder=3,marker='x')


# Add color bar
cb_ax_g = fig.add_axes([0.13, 0.10, 0.33, 0.01])
cb_g = fig.colorbar(tp_int, cax=cb_ax_g, orientation='horizontal')
cb_g.set_label('mm/year (%)', fontsize=16)
cb_g.ax.tick_params(labelsize=12)  # Set colorbar tick label font size to 14
# ------------------------ (h) MERRA2 Spatial Trend ------------------------
file_path = '/data/se_prev/prev_trend_and_p_1982_2020.nc'
nc_data = Dataset(file_path, 'r')

# Draw spatial distribution map of variable tp
var_name = 'trend'
variable = nc_data.variables[var_name][:]*100

print("Min:", np.min(variable))
print("Max:", np.max(variable))


# Get p_value data
p_value = nc_data.variables['p_value'][:]

# Get longitude and latitude
lon = nc_data.variables['longitude'][:]
lat = nc_data.variables['latitude'][:] 

# Create longitude and latitude grid
lon_grid, lat_grid = np.meshgrid(lon, lat)

# Draw map
ax_h = fig.add_subplot(4, 2, 8, projection=ccrs.PlateCarree())
# Create subplot in absolute coordinates of canvas
ax_h.set_extent([leftlon, rightlon, lowerlat, upperlat], crs=ccrs.PlateCarree())
# Define geographic coordinate label format
xticks = np.arange(leftlon, rightlon, 10)
yticks = np.arange(lowerlat, upperlat, 5)

# Remove first coordinate point
xticks = xticks[1:]  # Exclude first longitude
yticks = yticks[1:]  # Exclude first latitude

ax_h.set_xticks(xticks, crs=ccrs.PlateCarree())
ax_h.set_yticks(yticks, crs=ccrs.PlateCarree())
# Format longitude and latitude tick labels to standard format
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax_h.xaxis.set_major_formatter(lon_formatter)
ax_h.yaxis.set_major_formatter(lat_formatter)
ax_h.tick_params(axis='x', labelsize=16)  # X-axis font size
ax_h.tick_params(axis='y', labelsize=16)  # Y-axis font size
ax_h.set_title('Trend of precipitation variability in the SE',loc='center',fontsize=18)

# Read shp file
china = shpreader.Reader('/home/tangzhen/northwest/northwest_China.shp').geometries()
# Draw China national boundary with nine-dash line
ax_h.add_geometries(china, ccrs.PlateCarree(),facecolor='none', edgecolor='black', linewidth=1.5, zorder = 2)

# Set color bar breakpoints
boundaries1 = np.arange(-0.3, 0, 0.03)
boundaries2 = np.arange(0, 2.51 ,0.25)
boundaries = np.concatenate((boundaries1, boundaries2))
norm = BoundaryNorm(boundaries, ncolors=128)
# Draw spatial distribution map of tp-intensity
tp_int = ax_h.contourf(lon, lat, variable, cmap=cmaps.MPL_BrBG,extend='both', levels=boundaries,norm=norm, transform=ccrs.PlateCarree(),zorder=1)
#tp_int = f2_ax1.contourf(lon, lat,variable, cmap=cmaps.MPL_BrBG,extend='both',  transform=ccrs.PlateCarree())

# Mark grid points with p_value less than 0.05
significant_mask = p_value < 0.05
# Get longitude and latitude of significant grid points
significant_lons = lon_grid[significant_mask]
significant_lats = lat_grid[significant_mask]
# Mark significant grid points with small black dots on the map
ax_h.scatter(significant_lons, significant_lats, color='black', s=0.05, transform=ccrs.PlateCarree(), zorder=3,marker='x')


# Add color bar
cb_ax_h = fig.add_axes([0.56, 0.10, 0.33, 0.01])
cb_h = fig.colorbar(tp_int, cax=cb_ax_h, orientation='horizontal')
cb_h.set_label('mm/year (%)', fontsize=16)
cb_h.ax.tick_params(labelsize=12)  # Set colorbar tick label font size to 14

ax_a.text(0.00, 1.02, "(a)", transform=ax_a.transAxes, fontsize=20, ha='left', va='bottom', fontweight='bold')
ax_b.text(0.00, 1.02, "(b)", transform=ax_b.transAxes, fontsize=20, ha='left', va='bottom', fontweight='bold')
ax_c.text(0.00, 1.02, "(c)", transform=ax_c.transAxes, fontsize=20, ha='left', va='bottom', fontweight='bold')
ax_d.text(0.00, 1.02, "(d)", transform=ax_d.transAxes, fontsize=20, ha='left', va='bottom', fontweight='bold')
ax_e.text(0.00, 1.02, "(e)", transform=ax_e.transAxes, fontsize=20, ha='left', va='bottom', fontweight='bold')
ax_f.text(0.00, 1.02, "(f)", transform=ax_f.transAxes, fontsize=20, ha='left', va='bottom', fontweight='bold')
ax_g.text(0.00, 1.02, "(g)", transform=ax_g.transAxes, fontsize=20, ha='left', va='bottom', fontweight='bold')
ax_h.text(0.00, 1.02, "(h)", transform=ax_h.transAxes, fontsize=20, ha='left', va='bottom', fontweight='bold')


plt.savefig('', dpi=800, bbox_inches="tight")