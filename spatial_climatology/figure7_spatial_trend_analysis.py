import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature
from netCDF4 import Dataset
from matplotlib.colors import BoundaryNorm
import cartopy.mpl.ticker as cticker

leftlon, rightlon, lowerlat, upperlat = 70, 115, 30, 51
fig = plt.figure(figsize=(16, 12))
# Create adjusted blue-green-yellow-orange-red gradient, blue and green from light to dark
colors = [
    "#E5FFFF", "#B2FFFF", "#99FFFF", "#80FFFF", "#66CCFF", "#4D99FF", "#3366FF", "#1A33FF", "#0000FF", "#0000CC",  # Blue gradient (10 gradients) # Blue gradient (10 gradients, light to dark)
    "#CCFFCC", "#99FF99", "#66FF66", "#4DFF4D", "#33FF33", "#1AFF1A", "#00FF00", "#1AFF1A", "#FFFF66", "#FFFF33",  # Green gradient (10 gradients, light to deep)
    "#FFFA00", "#FFEA00", "#FFD700", "#FFCD00", "#FFB300", "#FF9900", "#FF8C00", "#FF7F00", "#FF6A00", "#FF5500",  # Yellow to orange (10 gradients)
    "#FF3300", "#FF1A00", "#FF0000", "#E60000", "#CC0000"  # Orange to red gradient (5 gradients)
]
# Create custom color gradient
n_bins = len(colors)  # Total number of color levels is 35 (including 5 new gradients from yellow to orange and 5 from orange to red)
cmap = LinearSegmentedColormap.from_list("custom_gradient", colors, N=n_bins)

boundaries1 = np.arange(0, 90, 3)
boundaries2 = np.arange(90, 101, 2)
boundaries = np.concatenate((boundaries1, boundaries2))
norm = BoundaryNorm(boundaries, ncolors=n_bins)
# ------------------------ (a) MERRA2 Spatial Trend ------------------------
file_path = '/data/nw_prev/nw_contributition1.nc'
nc_data = Dataset(file_path, 'r')

var_name = 'con'
variable = nc_data.variables[var_name][:]*100
print("Min:", np.min(variable))
print("Max:", np.max(variable))
lon = nc_data.variables['longitude'][:]
lat = nc_data.variables['latitude'][:] 
lon_grid, lat_grid = np.meshgrid(lon, lat)
ax_e = fig.add_subplot(2, 2, 1, projection=ccrs.PlateCarree())
ax_e.set_extent([leftlon, rightlon, lowerlat, upperlat], crs=ccrs.PlateCarree())

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
ax_e.set_title('Contribution of NW', loc='center', fontsize=20)

# Read shp file
china = shpreader.Reader('/home/tangzhen/northwest/northwest_China.shp').geometries()
# Draw China national boundary with nine-dash line
ax_e.add_geometries(china, ccrs.PlateCarree(),facecolor='none', edgecolor='black', linewidth=1.5, zorder = 2)

# Set color bar breakpoints
tp_int = ax_e.contourf(lon, lat, variable, cmap=cmap, extend='neither', levels=boundaries, norm=norm, transform=ccrs.PlateCarree(), zorder=1)
cb_ax_e = fig.add_axes([0.12, 0.52, 0.35, 0.015])
cb_e = fig.colorbar(tp_int, cax=cb_ax_e, orientation='horizontal')
cb_e.set_label('%', fontsize=16)
cb_e.ax.tick_params(labelsize=16)  # Set colorbar tick label font size to 12
ax_e.text(-0.05, 1.02, "(a)", transform=ax_e.transAxes, fontsize=18, ha='left', va='bottom', fontweight='bold')
# ------------------------ (b) MERRA2 Spatial Trend ------------------------
# Read nc file
file_path = '/data/ne_prev/ne_contributition1.nc'
nc_data = Dataset(file_path, 'r')
var_name = 'con'
variable = nc_data.variables[var_name][:]*100
print("Min:", np.min(variable))
print("Max:", np.max(variable))
lon = nc_data.variables['longitude'][:]
lat = nc_data.variables['latitude'][:] 
lon_grid, lat_grid = np.meshgrid(lon, lat)

# Draw map
ax_f = fig.add_subplot(2, 2, 2, projection=ccrs.PlateCarree())
ax_f.set_extent([leftlon, rightlon, lowerlat, upperlat], crs=ccrs.PlateCarree())
xticks = np.arange(leftlon, rightlon, 10)
yticks = np.arange(lowerlat, upperlat, 5)
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
ax_f.set_title('Contribution of NE', loc='center', fontsize=20)

# Read shp file
china = shpreader.Reader('/home/tangzhen/northwest/northwest_China.shp').geometries()
ax_f.add_geometries(china, ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=1.5, zorder=2)

# Set color bar breakpoints
tp_int = ax_f.contourf(lon, lat, variable, cmap=cmap, extend='neither', levels=boundaries, norm=norm, transform=ccrs.PlateCarree(), zorder=1)

cb_ax_f = fig.add_axes([0.55, 0.52, 0.35, 0.015])
cb_f = fig.colorbar(tp_int, cax=cb_ax_f, orientation='horizontal')
cb_f.set_label('%', fontsize=16)
cb_f.ax.tick_params(labelsize=16)  # Set colorbar tick label font size to 12
ax_f.text(-0.05, 1.02, "(b)", transform=ax_f.transAxes, fontsize=18, ha='left', va='bottom', fontweight='bold')
# ------------------------ (c) MERRA2 Spatial Trend ------------------------
# Read nc file
file_path = '/data/sw_prev/sw_contributition.nc'
nc_data = Dataset(file_path, 'r')

# Draw spatial distribution map of variable tp
var_name = 'con'
variable = nc_data.variables[var_name][:]*100

print("Min:", np.min(variable))
print("Max:", np.max(variable))


# Get longitude and latitude
lon = nc_data.variables['longitude'][:]
lat = nc_data.variables['latitude'][:] 

# Create longitude and latitude grid
lon_grid, lat_grid = np.meshgrid(lon, lat)
ax_g = fig.add_subplot(2, 2, 3, projection=ccrs.PlateCarree())
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
ax_g.set_title('Contribution of SW', loc='center', fontsize=20)

# Read shp file
china = shpreader.Reader('/home/tangzhen/northwest/northwest_China.shp').geometries()
# Draw China national boundary with nine-dash line
ax_g.add_geometries(china, ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=1.5, zorder=2)

# Set color bar breakpoints
tp_int = ax_g.contourf(lon, lat, variable, cmap=cmap, extend='neither', levels=boundaries, norm=norm, transform=ccrs.PlateCarree(), zorder=1)

cb_ax_g = fig.add_axes([0.12, 0.10, 0.35, 0.015])
cb_g = fig.colorbar(tp_int, cax=cb_ax_g, orientation='horizontal')
cb_g.set_label('%', fontsize=16)
cb_g.ax.tick_params(labelsize=16)  # Set colorbar tick label font size to 12
ax_g.text(-0.05, 1.02, "(c)", transform=ax_g.transAxes, fontsize=18, ha='left', va='bottom', fontweight='bold')
# ------------------------ (d) MERRA2 Spatial Trend ------------------------
file_path = '/data/se_prev/se_contributition.nc'
nc_data = Dataset(file_path, 'r')

# Draw spatial distribution map of variable tp
var_name = 'con'
variable = nc_data.variables[var_name][:]*100

print("Min:", np.min(variable))
print("Max:", np.max(variable))


# Get longitude and latitude
lon = nc_data.variables['longitude'][:]
lat = nc_data.variables['latitude'][:] 

# Create longitude and latitude grid
lon_grid, lat_grid = np.meshgrid(lon, lat)
ax_h = fig.add_subplot(2, 2, 4, projection=ccrs.PlateCarree())
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
ax_h.set_title('Contribution of SE', loc='center', fontsize=20)

# Read shp file
china = shpreader.Reader('/home/tangzhen/northwest/northwest_China.shp').geometries()
# Draw China national boundary with nine-dash line
ax_h.add_geometries(china, ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=1.5, zorder=2)


# Draw spatial distribution map of tp-intensity
tp_int = ax_h.contourf(lon, lat, variable, cmap=cmap, extend='neither', levels=boundaries, norm=norm, transform=ccrs.PlateCarree(), zorder=1)

# Add color bar
cb_ax_h = fig.add_axes([0.55, 0.10, 0.35, 0.015])
cb_h = fig.colorbar(tp_int, cax=cb_ax_h, orientation='horizontal')
cb_h.set_label('%', fontsize=16)
cb_h.ax.tick_params(labelsize=16)  # Set colorbar tick label font size to 12
ax_h.text(-0.05, 1.02, "(d)", transform=ax_h.transAxes, fontsize=18, ha='left', va='bottom', fontweight='bold')

plt.savefig('', dpi=800, bbox_inches="tight")