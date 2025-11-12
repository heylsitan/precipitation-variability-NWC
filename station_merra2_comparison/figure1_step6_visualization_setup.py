import os
import pandas as pd 
import numpy as np
import geopandas as gpd
import xarray as xr
from scipy.stats import pearsonr
from scipy.stats import linregress
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy.spatial import Voronoi
from shapely.geometry import Point, Polygon
from sklearn.neighbors import NearestNeighbors

def voronoi_finite_polygons_2d(vor, radius=1e6):
    """Process infinite Voronoi polygon regions"""
    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]
        if all(v >= 0 for v in vertices):
            new_regions.append(vertices)
            continue

        ridges = all_ridges.get(p1, [])
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                continue

            t = vor.points[p2] - vor.points[p1]
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_vertices.append(far_point.tolist())
            new_region.append(len(new_vertices) - 1)

        new_regions.append(new_region)

    return new_regions, np.array(new_vertices)

# Read station longitude and latitude
sll_file_path = '/data/pre_data/station_lon_lat.xlsx'
station_lat_lon = pd.read_excel(sll_file_path)
station_numbers = station_lat_lon['stationnumber']
lats = station_lat_lon['latitude']
lons = station_lat_lon['longitude']

# Read Northwest China boundary
shp_file = '/home/northwest/northwest_China.shp'
nwc_shp = gpd.read_file(shp_file)

# Filter stations within Northwest China
nwc_stations = []
valid_indices = []  # Record indices of valid stations
for i in range(len(station_numbers)):
    point = Point(lons[i], lats[i])
    if nwc_shp.contains(point).any():
        nwc_stations.append(station_numbers[i])
        valid_indices.append(i)

# Extract longitude and latitude of valid stations
valid_lons = [lons[i] for i in valid_indices]
valid_lats = [lats[i] for i in valid_indices]
# ------------------------- Voronoi Polygon Weight Calculation -------------------------
# Generate Voronoi polygons
points = np.array([[lon, lat] for lon, lat in zip(valid_lons, valid_lats)])
vor = Voronoi(points)
regions, vertices = voronoi_finite_polygons_2d(vor)

# Build GeoDataFrame
voronoi_polys = [Polygon(vertices[region]) for region in regions]
voronoi_gdf = gpd.GeoDataFrame(
    geometry=voronoi_polys,
    crs=nwc_shp.crs  # Ensure coordinate system consistency
)

# Clip to Northwest China
nwc_voronoi = gpd.overlay(voronoi_gdf, nwc_shp, how='intersection')

# Generate station GeoDataFrame
stations_gdf = gpd.GeoDataFrame(
    geometry=gpd.points_from_xy(valid_lons, valid_lats),
    crs=nwc_shp.crs
)
stations_gdf['station_id'] = nwc_stations

# Calculate Voronoi polygon centroids for spatial matching
nwc_voronoi['centroid'] = nwc_voronoi.geometry.centroid

# Perform nearest neighbor spatial join
nwc_voronoi = gpd.sjoin_nearest(
    nwc_voronoi.to_crs(stations_gdf.crs),
    stations_gdf[['geometry', 'station_id']],
    how='left',
    distance_col="distance"
)

# Clean redundant columns and calculate weights
nwc_voronoi = nwc_voronoi[['geometry', 'station_id']]
total_area = nwc_voronoi.geometry.area.sum()
nwc_voronoi['area_weight'] = nwc_voronoi.geometry.area / total_area
# ------------------------- Precipitation Data Processing -------------------------
# Read precipitation data
pre_file_path = '/data/pre_data/station_pre.xlsx'
df = pd.read_excel(pre_file_path)
df['Times'] = pd.to_datetime(df.iloc[:,0])
df['year'] = df['Times'].dt.year
df['month'] = df['Times'].dt.month

# Time range settings
station_years = np.arange(1982, 2019)
summer_months = [6, 7, 8]  # Define summer months

# Initialize storage structure
summer_precip = {year: {} for year in station_years}

# Calculate summer total precipitation for each station
for year in station_years:
    print(f"Processing {year}...")
    for station_id in nwc_stations:
        if station_id not in df.columns:
            continue
            
        # Extract summer data
        summer_data = df[
            (df['year'] == year) & 
            (df['month'].isin(summer_months))
        ][station_id]
        
        # Summer total precipitation (ignore missing values)
        total_precip = summer_data.std(ddof=0, skipna=True)
        summer_precip[year][station_id] = total_precip


# ------------------------- Missing Value Interpolation -------------------------
# Fill missing values using nearest neighbor station values
coords = np.array([[lon, lat] for lon, lat in zip(valid_lons, valid_lats)])
nbrs = NearestNeighbors(n_neighbors=1).fit(coords)

for year in station_years:
    for idx, station_id in enumerate(nwc_stations):
        if np.isnan(summer_precip[year].get(station_id, np.nan)):
            # Find nearest station
            _, nearest_idx = nbrs.kneighbors([coords[idx]])
            nearest_id = nwc_stations[nearest_idx[0][0]]
            summer_precip[year][station_id] = summer_precip[year][nearest_id]

# ------------------------- Weighted Average Calculation -------------------------
mean_value = []

for year in station_years:
    print(year)
    weighted_sum = 0.0
    for _, row in nwc_voronoi.iterrows():
        station_id = row['station_id']
        precip = summer_precip[year].get(station_id, np.nan)
        if not np.isnan(precip):
            weighted_sum += precip * row['area_weight']

    mean_value.append(weighted_sum)

print(mean_value)

# Save as Excel
df_mean = pd.DataFrame({
    'Year': station_years,
    'Prev': mean_value
})
save_path = '/data/station_prev.xlsx'
df_mean.to_excel(save_path, index=False)
print(f"Saved successfully：{save_path}")