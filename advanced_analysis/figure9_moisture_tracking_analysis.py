import xarray as xr
import os
import numpy as np
import cmaps
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from netCDF4 import Dataset
from matplotlib.colors import BoundaryNorm
import cartopy.mpl.ticker as cticker
import matplotlib.gridspec as gridspec

leftlon, rightlon, lowerlat, upperlat = 70, 115, 30, 50
# === Set up figure ===
fig = plt.figure(figsize=(24, 24))
gs = gridspec.GridSpec(2, 2, height_ratios=[2, 1])  # Top plots are larger (2), bottom ones smaller (1)

# ------------------------ (a)------------------------
file_path = ''
nc_data = Dataset(file_path, 'r')
var_name = ''
variable = nc_data.variables[var_name][:]
lon = nc_data.variables['longitude'][:]
lat = nc_data.variables['latitude'][:]

print("Min:", np.min(variable))
print("Max:", np.max(variable))

# === Set up image and projection ===
ax_a = fig.add_subplot(gs[0, 0], projection=ccrs.PlateCarree())  # Top-left plot
# Set display area range
ax_a.set_extent([-30, 180, -50, 90], crs=ccrs.PlateCarree())  # Longitude range [60, 135], latitude range [15, 55]
ax_a.coastlines()
#ax_a.gridlines(draw_labels=True)
xticks = np.arange(0, 180, 30)  # Longitude ticks
yticks = np.arange(-40, 90, 20)  # Latitude ticks

# Set x and y axis ticks
ax_a.set_xticks(xticks)
ax_a.set_yticks(yticks)

# Set x and y axis labels
ax_a.set_xticklabels([f'{x}°E' for x in xticks], fontsize=16)  # Manually set longitude labels
yticklabels = [f'{abs(y)}°{"N" if y >= 0 else "S"}' for y in yticks]
ax_a.set_yticklabels(yticklabels, fontsize=16)  # Manually set latitude labels, dynamically add N/S

ax_a.grid(True, color='gray', linestyle='-', linewidth=0.5)

# === Set color bar breakpoints ===
boundaries1 = np.arange(-1.0, 0.0, 1.0)
boundaries2 = np.arange(0.1, 10.0 ,2.0)
boundaries3 = np.arange(10.0,35.0, 5.0)
boundaries = np.concatenate((boundaries1, boundaries2, boundaries3))
norm = BoundaryNorm(boundaries, ncolors=17)

tp_plot = ax_a.contourf(lon, lat,variable, cmap=cmaps.precip3_16lev,extend='both', levels=boundaries,norm=norm, transform=ccrs.PlateCarree())
#tp_plot = ax_a.contourf(lon, lat,variable, cmap=cmaps.MPL_BrBG,extend='both', transform=ccrs.PlateCarree())
NW_shp = '/data/shp/NW_region.shp'  # If you have a more complete China boundary shp, that would be better
NW_geom = shpreader.Reader(NW_shp).geometries()
ax_a.add_geometries(NW_geom, ccrs.PlateCarree(),
                  facecolor='none', edgecolor='#CA0E12', linewidth=1.5, zorder=3)

NE_shp = '/data/shp/NE_region.shp'  # If you have a more complete China boundary shp, that would be better
NE_geom = shpreader.Reader(NE_shp).geometries()
ax_a.add_geometries(NE_geom, ccrs.PlateCarree(),
                  facecolor='none', edgecolor='#25377F', linewidth=1.5, zorder=3)

SW_shp = '/data/shp/SW_region.shp'  # If you have a more complete China boundary shp, that would be better
SW_geom = shpreader.Reader(SW_shp).geometries()
ax_a.add_geometries(SW_geom, ccrs.PlateCarree(),
                  facecolor='none', edgecolor='#2AA7DE', linewidth=1.5, zorder=3)

SE_shp = '/data/shp/SE_region.shp'  # If you have a more complete China boundary shp, that would be better
SE_geom = shpreader.Reader(SE_shp).geometries()
ax_a.add_geometries(SE_geom, ccrs.PlateCarree(),
                  facecolor='none', edgecolor='#F6BD21', linewidth=1.5, zorder=3)

ax_a.text(50, 60, 'NW', transform=ccrs.PlateCarree(),
        fontsize=16, fontweight='bold', color='#CA0E12', ha='center')

ax_a.text(110, 60, 'NE', transform=ccrs.PlateCarree(),
        fontsize=16, fontweight='bold', color='#25377F', ha='center')

ax_a.text(62, 15, 'SW', transform=ccrs.PlateCarree(),
        fontsize=16, fontweight='bold', color='#2AA7DE', ha='center')

ax_a.text(125, 15, 'SE', transform=ccrs.PlateCarree(),
        fontsize=16, fontweight='bold', color='#F6BD21', ha='center')
ax_a.text(90, 38, 'NWC', transform=ccrs.PlateCarree(),
        fontsize=16, fontweight='bold', color='black', ha='center')

# Read shp file
china = shpreader.Reader('/home/tangzhen/northwest/northwest_China.shp').geometries()
# Draw China border with nine-dash line
ax_a.add_geometries(china, ccrs.PlateCarree(),facecolor='none', edgecolor='black', linewidth=1.5,zorder=5)

# === Add color bar ===
# === Add color bar ===
cb_ax_a = fig.add_axes([0.12, 0.46, 0.35, 0.01])
cb_a = fig.colorbar(tp_plot, cax=cb_ax_a, orientation='horizontal')
#cbar.set_ticks(np.arange(0.0, 7.2, 0.8)) # Only show major ticks
cb_a.ax.tick_params(which='minor', length=0)  # Remove minor tick lines
cb_a.set_label('mm', fontsize=18)
cb_a.ax.tick_params(labelsize=16)  # Set colorbar tick label font size to 14
ax_a.set_title("Difference in moisture source", fontsize=20, loc='center')
ax_a.text(0.00, 1.02, "(a)", transform=ax_a.transAxes, fontsize=20, ha='left', va='bottom', fontweight='bold')
# ------------------------ (b)------------------------
# Read nc file
file_path = '/data/tangzhen/difference.nc'
nc_data = Dataset(file_path, 'r')

# Draw spatial distribution of variable tp
var_name = 'tp'
variable = nc_data.variables[var_name][:]*1000
print(f"Min: {np.min(variable):.8f}")  # Keep 2 decimal places
print(f"Max: {np.max(variable):.8f}")  # Keep 2 decimal places

# Get longitude and latitude
lon = nc_data.variables['lon'][:]
lat = nc_data.variables['lat'][:]

ax_b = fig.add_subplot(gs[0, 1], projection=ccrs.PlateCarree())  # Top-right plot
ax_b.set_extent([leftlon, rightlon, lowerlat, upperlat], crs=ccrs.PlateCarree())
# Define geographic coordinate label format
xticks = np.arange(leftlon, rightlon, 10)
yticks = np.arange(lowerlat, upperlat, 5)

# Remove first coordinate point
xticks = xticks[1:]  # Exclude first longitude
yticks = yticks[1:]  # Exclude first latitude

ax_b.set_xticks(xticks, crs=ccrs.PlateCarree())
ax_b.set_yticks(yticks, crs=ccrs.PlateCarree())
# Format longitude and latitude tick labels in standard format
lon_formatter = cticker.LongitudeFormatter()
lat_formatter = cticker.LatitudeFormatter()
ax_b.xaxis.set_major_formatter(lon_formatter)
ax_b.yaxis.set_major_formatter(lat_formatter)
ax_b.tick_params(axis='x', labelsize=16)  # x-axis font size
ax_b.tick_params(axis='y', labelsize=16)  # y-axis font size
ax_b.set_title('Difference in precipitation variability', loc='center', fontsize=20)

china = shpreader.Reader('/home/tangzhen/northwest/northwest_China.shp').geometries()
ax_b.add_geometries(china, ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=1.5, zorder=2)

boundaries = np.arange(0.0, 4.1, 0.1)
norm = BoundaryNorm(boundaries, ncolors=128)
tp_int = ax_b.contourf(lon, lat, variable, cmap=cmaps.MPL_YlGnBu, extend='both', levels=boundaries, norm=norm, transform=ccrs.PlateCarree(), zorder=1)

cb_ax_b = fig.add_axes([0.55, 0.46, 0.35, 0.01])
cb_b = fig.colorbar(tp_int, cax=cb_ax_b, orientation='horizontal')
cb_b.set_label('mm', fontsize=18)
cb_b.ax.tick_params(labelsize=16)  # Set colorbar tick label font size to 14

ax_b.text(0.00, 1.02, "(b)", transform=ax_b.transAxes, fontsize=20, ha='left', va='bottom', fontweight='bold')
# ------------------------ (c)------------------------
# === Dry and wet years ===
dry_years = np.array([1983, 1985, 1986, 1994, 1995])
wet_years = np.array([2005, 2007, 2012, 2013, 2016, 2017, 2018])

# === Regional paths ===
paths = {
    'NW': '/data/nw_prev',
    'NE': '/data/ne_prev',
    'SW': '/data/sw_prev',
    'SE': '/data/se_prev',
    'NWC': '/data/internal_prev'
}

# === Keep consistent region order and colors ===
regions = ['NWC', 'NW', 'NE', 'SW', 'SE']
colors = ['dimgrey','#CA0E12', '#25377F', '#2AA7DE', '#F6BD21']

# === Function to read precipitation data ===
def read_total_tp(path, prefix, years):
    yearly_totals = []
    for year in years:
        file_path = os.path.join(path, f'{prefix}_{year}.nc')
        ds = xr.open_dataset(file_path)
        data = ds['tp']
        if prefix == 'internal_prev':
            data = data * 1000  # Unit conversion
        lat_name = 'lat' if 'lat' in data.dims else 'latitude'
        lon_name = 'lon' if 'lon' in data.dims else 'longitude'
        weights = np.cos(np.deg2rad(data[lat_name]))
        weighted_avg = data.weighted(weights).mean(dim=[lat_name, lon_name])
        yearly_totals.append(float(weighted_avg.sum().values))
    return np.mean(yearly_totals)

# === Get dry and wet year data ===
dry_values = [
    read_total_tp(paths[r], 'prev', dry_years) if r != 'NWC' else read_total_tp(paths[r], 'internal_prev', dry_years)
    for r in regions
]
wet_values = [
    read_total_tp(paths[r], 'prev', wet_years) if r != 'NWC' else read_total_tp(paths[r], 'internal_prev', wet_years)
    for r in regions
]

# Calculate residuals
total_path = '/data/tangzhen/merra2_prev_nwc'
# Get paths for all .nc files
dry_file_paths = [os.path.join(total_path, f'prev_{year}.nc') for year in dry_years]
# Calculate total precipitation variability for low years
dry_prev_mean = []
for file_path in dry_file_paths:
    ds = xr.open_dataset(file_path)
    data = ds['tp']*1000
    weights = np.cos(np.deg2rad(data.lat))
    weitght_avg = data.weighted(weights).mean(dim=["lat","lon"])
    dry_prev_mean.append(weitght_avg)

dry_mean_value = np.mean(dry_prev_mean)

wet_file_paths = [os.path.join(total_path, f'prev_{year}.nc') for year in wet_years]
# Calculate total precipitation variability for high years
wet_prev_mean = []
for file_path in wet_file_paths:
    ds = xr.open_dataset(file_path)
    data = ds['tp']*1000
    weights = np.cos(np.deg2rad(data.lat))
    weitght_avg = data.weighted(weights).mean(dim=["lat","lon"])
    wet_prev_mean.append(weitght_avg)

wet_mean_value = np.mean(wet_prev_mean)

# Calculate residuals for dry and wet years
total_dry = sum(dry_values)  # Total precipitation for dry years
total_wet = sum(wet_values)  # Total precipitation for wet years

# Calculate residuals for dry and wet years
residual_dry = dry_mean_value - total_dry  # Dry year residual = total precipitation - dry year precipitation
residual_wet = wet_mean_value - total_wet  # Wet year residual = total precipitation - wet year precipitation


# === Calculate increase amounts and percentages ===
increase_amounts = [wet - dry for wet, dry in zip(wet_values, dry_values)]
res_increase = residual_wet - residual_dry
total_increase = sum(increase_amounts) + res_increase
increase_percentages = [(inc / total_increase) * 100 for inc in increase_amounts]
increase_percentages.append((res_increase/total_increase)*100)

# Calculate precipitation difference between wet and dry years for each region
differences = [wet - dry for wet, dry in zip(wet_values, dry_values)]

# Output precipitation difference for each region
for region, diff in zip(regions, differences):
    print(f"Difference for region {region}: {diff:.2f} mm")
# === Figure 1: Stacked bar chart for dry and wet years ===
ax_c = fig.add_subplot(gs[1, 0])  # Bottom-left plot

x = [-0.1, 0.1]
width = 0.04
bottom_dry = 0
bottom_wet = 0

for i in range(len(regions)):
    ax_c.bar(x[0], dry_values[i], width, bottom=bottom_dry, color=colors[i], edgecolor='black')
    ax_c.bar(x[1], wet_values[i], width, bottom=bottom_wet, color=colors[i], edgecolor='black')
    bottom_dry += dry_values[i]
    bottom_wet += wet_values[i]

# Add residual term for dry years
ax_c.bar(x[0], residual_dry, width, bottom=bottom_dry, color='mediumseagreen', edgecolor='black', label='RES')
bottom_dry += residual_dry  # Update dry year bottom position

# Add residual term for wet years
ax_c.bar(x[1], residual_wet, width, bottom=bottom_wet, color='mediumseagreen', edgecolor='black', label='RES')
bottom_wet += residual_wet  # Update wet year bottom position


ax_c.set_xticks(x)
ax_c.set_xticklabels(['Low', 'High'], fontsize=16)
ax_c.set_ylabel('Precipitation Variability (mm)', fontsize=16)
ax_c.set_xlim(-0.2, 0.2)
ax_c.set_ylim(0, max(bottom_dry, bottom_wet) * 1.2)
ax_c.tick_params(axis='y', labelsize=16)

legend_patches = [plt.Rectangle((0, 0), 1, 1, facecolor=colors[i], edgecolor='black') for i in range(len(regions))]
legend_patches.append(plt.Rectangle((0, 0), 1, 1, facecolor='mediumseagreen', edgecolor='black'))  # Add legend for residual
legend = ax_c.legend(
    legend_patches,
    regions + ['RES'],
    loc='upper left',
    bbox_to_anchor=(0.00, 0.98),
    fontsize=16,
    ncol=1,
    frameon=False
)
# Set legend border color to gray
legend.get_frame().set_edgecolor('gray')
legend.get_frame().set_linewidth(1.2)  # Optional: set border width
ax_c.set_title('Absolute contribution',loc='center',fontsize=20)
ax_c.text(0.01, 1.02, "(c)", transform=ax_c.transAxes, fontsize=20, ha='left', va='bottom', fontweight='bold')

# === Figure 2: Pie chart for regional precipitation increase proportion ===
ax_d = fig.add_subplot(gs[1, 1])  # Bottom-right plot
regions = regions + ['RES']

ax_d.pie(
    increase_percentages,
    labels=regions,
    colors=colors + ['mediumseagreen'],
    autopct='%1.1f%%',
    startangle=90,
    wedgeprops={'edgecolor': 'white'},
    textprops={'fontsize': 24}  # Set font size
)
ax_d.set_title('Relative contribution', fontsize=20, y=1.02)

ax_d.text(0.00, 1.02, "(d)", transform=ax_d.transAxes, fontsize=20, ha='left', va='bottom', fontweight='bold')
# ------------------------ (b)------------------------
plt.subplots_adjust(hspace=0.005)  # Default is 0.5, reduce it to make top-bottom closer

plt.savefig('', dpi=600, bbox_inches="tight")