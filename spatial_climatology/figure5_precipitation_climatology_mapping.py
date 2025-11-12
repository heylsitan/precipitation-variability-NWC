import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy.stats import linregress

# === Set up canvas ===
fig = plt.figure(figsize=(24, 24))

#=== (a)NW ====
folder_path = '/data/nw_prev'
file_paths = [os.path.join(folder_path, f'prev_{year}.nc') for year in range(1982, 2021)]
pre_int_mean = []
for file_path in file_paths:
    ds = xr.open_dataset(file_path)
    data = ds['tp']
    weights = np.cos(np.deg2rad(data.latitude))
    weitght_avg = data.weighted(weights).mean(dim=["latitude","longitude"])
    pre_int_mean.append(weitght_avg)

mean_value = np.mean(pre_int_mean)
years = np.arange(1982,2021)
slop, intercept, r_value, p_value, std_err = linregress(years, pre_int_mean)
trend_line = slop*years + intercept

ax_a = fig.add_subplot(3, 2, 1)
ax_a.plot(years, pre_int_mean, marker='s', color='#CA0E12', markersize=10, linewidth=2, linestyle='--')
ax_a.plot(years, trend_line, linestyle='-', color='black', linewidth=2)
ax_a.xaxis.set_minor_locator(MultipleLocator(1))
ax_a.xaxis.set_major_locator(MultipleLocator(5))
ax_a.tick_params(axis='both',
                which='minor',
                width=1,
                length=4,
                labelsize=16)
ax_a.tick_params(axis='both',
                which='major',
                width=1,
                length=6,
                labelsize=16)
if p_value < 0.05:
    ax_a.annotate(f'trend={slop:.4f}mm/year (p<0.05)',
                 xy=(years[5], trend_line[5]),
                 xytext=(-15, 210),
                 textcoords='offset points',
                 color='#CA0E12',
                 fontsize=18)
else:
    ax_a.annotate(f'trend={slop:.4f}mm/year(p={p_value:.2f})',
                 xy=(years[0], trend_line[10]),
                 xytext=(-10, 200),
                 textcoords='offset points',
                 color='#CA0E12',
                 fontsize=20)

ax_a.set_title('Precipitation variability from the NW',loc='center',fontsize=20)
ax_a.set_ylabel('mm',fontsize=18)
ax_a.set_xlabel('Years',fontsize=18)
ax_a.text(0.00, 1.02, "(a)", transform=plt.gca().transAxes, fontsize=20, ha='left', va='bottom', fontweight='bold')

# ===(b)NE ===
folder_path = '/data/ne_prev'
file_paths = [os.path.join(folder_path, f'prev_{year}.nc') for year in range(1982, 2021)]
pre_int_mean = []
for file_path in file_paths:
    ds = xr.open_dataset(file_path)
    data = ds['tp']
    weights = np.cos(np.deg2rad(data.latitude))
    weitght_avg = data.weighted(weights).mean(dim=["latitude","longitude"])
    pre_int_mean.append(weitght_avg)

mean_value = np.mean(pre_int_mean)
years = np.arange(1982,2021)
slop, intercept, r_value, p_value, std_err = linregress(years, pre_int_mean)
trend_line = slop*years + intercept

ax_b = fig.add_subplot(3, 2, 2)
ax_b.plot(years, pre_int_mean, marker='s', color='#25377F', markersize=10, linewidth=2, linestyle='--')
ax_b.plot(years, trend_line, linestyle='-', color='black', linewidth=2)
ax_b.xaxis.set_minor_locator(MultipleLocator(1))
ax_b.xaxis.set_major_locator(MultipleLocator(5))
ax_b.tick_params(axis='both',
                 which='minor',
                 width=1,
                 length=4,
                 labelsize=16)
ax_b.tick_params(axis='both',
                 which='major',
                 width=1,
                 length=6,
                 labelsize=16)

if p_value < 0.05:
    ax_b.annotate(f'trend={slop:.4f}mm/year (p<0.05)',
                  xy=(years[5], trend_line[5]),
                  xytext=(-80, 260),
                  textcoords='offset points',
                  color='#25377F',
                  fontsize=18)
else:
    ax_b.annotate(f'trend={slop:.4f}mm/year(p={p_value:.2f})',
                  xy=(years[5], trend_line[5]),
                  xytext=(-80, 200),
                  textcoords='offset points',
                  color='#25377F',
                  fontsize=20)

ax_b.set_title('Precipitation variability from the NE', loc='center', fontsize=20)
ax_b.set_ylabel('mm', fontsize=18)
ax_b.set_xlabel('Years', fontsize=18)
ax_b.text(0.00, 1.02, "(b)", transform=ax_b.transAxes, fontsize=20, ha='left', va='bottom', fontweight='bold')

# ===(c)SW ===
folder_path = '/data/sw_prev'
file_paths = [os.path.join(folder_path, f'prev_{year}.nc') for year in range(1982, 2021)]
pre_int_mean = []
for file_path in file_paths:
    ds = xr.open_dataset(file_path)
    data = ds['tp']
    weights = np.cos(np.deg2rad(data.latitude))
    weitght_avg = data.weighted(weights).mean(dim=["latitude","longitude"])
    pre_int_mean.append(weitght_avg)

mean_value = np.mean(pre_int_mean)
years = np.arange(1982,2021)
slop, intercept, r_value, p_value, std_err = linregress(years, pre_int_mean)
trend_line = slop*years + intercept

ax_c = fig.add_subplot(3,2,3)
ax_c.plot(years, pre_int_mean, marker='s', color='#2AA7DE', markersize=10, linewidth=2, linestyle='--')
ax_c.plot(years, trend_line, linestyle='-', color='black', linewidth=2)
ax_c.xaxis.set_minor_locator(MultipleLocator(1))
ax_c.xaxis.set_major_locator(MultipleLocator(5))
ax_c.tick_params(axis='both',
                 which='minor',
                 width=1,
                 length=4,
                 labelsize=16)
ax_c.tick_params(axis='both',
                 which='major',
                 width=1,
                 length=6,
                 labelsize=16)


if p_value < 0.05:
    ax_c.annotate(f'trend={slop:.4f}mm/year (p<0.05)',
                  xy=(years[5], trend_line[5]),
                  xytext=(-15, 220),
                  textcoords='offset points',
                  color='#2AA7DE',
                  fontsize=20)
else:
    ax_c.annotate(f'trend={slop:.4f}mm/year(p={p_value:.2f})',
                  xy=(years[3], trend_line[3]),
                  xytext=(-15, 120),
                  textcoords='offset points',
                  color='#2AA7DE',
                  fontsize=18)

ax_c.set_title('Precipitation variability from the SW', loc='center', fontsize=20)
ax_c.set_ylabel('mm', fontsize=18)
ax_c.set_xlabel('Years', fontsize=18)
ax_c.text(0.00, 1.02, "(c)", transform=ax_c.transAxes, fontsize=20, ha='left', va='bottom', fontweight='bold')

# ==(d)SE ===
folder_path = '/data/se_prev'
file_paths = [os.path.join(folder_path, f'prev_{year}.nc') for year in range(1982, 2021)]
pre_int_mean = []
for file_path in file_paths:
    ds = xr.open_dataset(file_path)
    data = ds['tp']
    weights = np.cos(np.deg2rad(data.latitude))
    weitght_avg = data.weighted(weights).mean(dim=["latitude","longitude"])
    pre_int_mean.append(weitght_avg)

mean_value = np.mean(pre_int_mean)
years = np.arange(1982,2021)
slop, intercept, r_value, p_value, std_err = linregress(years, pre_int_mean)
trend_line = slop*years + intercept

ax_d = fig.add_subplot(3, 2, 4)
ax_d.plot(years, pre_int_mean, marker='s', color='#F6BD21', markersize=10, linewidth=2, linestyle='--')
ax_d.plot(years, trend_line, linestyle='-', color='black', linewidth=2)
ax_d.xaxis.set_minor_locator(MultipleLocator(1))
ax_d.xaxis.set_major_locator(MultipleLocator(5))
ax_d.tick_params(axis='both',
                 which='minor',
                 width=1,
                 length=4,
                 labelsize=16)
# Set major tick style
ax_d.tick_params(axis='both',
                 which='major',
                 width=1,
                 length=6,
                 labelsize=16)
if p_value < 0.05:
    ax_d.annotate(f'trend={slop:.4f}mm/year (p<0.05)',
                  xy=(years[5], trend_line[5]),
                  xytext=(-15, 200),
                  textcoords='offset points',
                  color='#F6BD21',
                  fontsize=20)
else:
    ax_d.annotate(f'trend={slop:.4f}mm/year(p={p_value:.2f})',
                  xy=(years[5], trend_line[5]),
                  xytext=(0, 200),
                  textcoords='offset points',
                  color='#F6BD21',
                  fontsize=18)

ax_d.set_title('Precipitation variability from the SE', loc='center', fontsize=20)
ax_d.set_ylabel('mm', fontsize=18)
ax_d.set_xlabel('Years', fontsize=18)
ax_d.text(0.00, 1.02, "(d)", transform=ax_d.transAxes, fontsize=20, ha='left', va='bottom', fontweight='bold')

# ===(e)total amount ===
years = np.arange(1982, 2021)
paths = {
    'NW': '/data/nw_prev',
    'NE': '/data/ne_prev',
    'SW': '/data/sw_prev',
    'SE': '/data/se_prev'
}

def read_total_tp(path, prefix):
    yearly_totals = []
    for year in years:
        file_path = os.path.join(path, f'{prefix}_{year}.nc')
        ds = xr.open_dataset(file_path)
        data = ds['tp']
        lat_name = 'lat' if 'lat' in data.dims else 'latitude'
        lon_name = 'lon' if 'lon' in data.dims else 'longitude'
        weights = np.cos(np.deg2rad(data[lat_name]))
        weighted_avg = data.weighted(weights).mean(dim=[lat_name, lon_name])
        total_in_year = float(weighted_avg.sum().values)
        yearly_totals.append(total_in_year)
    return np.mean(yearly_totals)

total_amounts = {
    region: read_total_tp(path, "prev")
    for region, path in paths.items()
}

total_sum_all = sum(total_amounts.values())
contributions = {k: v / total_sum_all for k, v in total_amounts.items()}

labels = list(contributions.keys())
sizes = [contributions[k]*100 for k in labels]  # Percentage
colors = ['#CA0E12', '#25377F', '#2AA7DE', '#F6BD21']

ax_e = fig.add_subplot(3,2,5)
bars = ax_e.bar(labels, sizes, width=0.4, color=colors, edgecolor='black')
for bar, size in zip(bars, sizes):
    ax_e.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5, f'{size:.1f}%', 
              ha='center', va='bottom', fontsize=20)

# Beautify chart
ax_e.set_ylim(0, max(sizes)*1.15)
ax_e.set_ylabel('Percentage (%)', fontsize=18)
ax_e.set_title('Relative contribution to total amuont', fontsize=20)
ax_e.set_xticklabels(labels, fontsize=16)
ax_e.set_yticklabels(ax_e.get_yticks(), fontsize=16)
ax_e.text(0.00, 1.02, "(e)", transform=ax_e.transAxes, fontsize=20, ha='left', va='bottom', fontweight='bold')

# ===Trend ===
# ===== Year range =====
years = np.arange(1982, 2021)

# ===== Path dictionary (unified management) =====
paths = {
    'NW': '/data/nw_prev',
    'NE': '/data/prev',
    'SW': '/data/prev',
    'SE': '/data/se_prev'
}

# ===== Reading function =====
def read_weighted_tp(path, prefix):
    tp_list = []
    for year in years:
        file_path = os.path.join(path, f'{prefix}_{year}.nc')
        ds = xr.open_dataset(file_path)
        data = ds['tp']
        lat_name = 'lat' if 'lat' in data.dims else 'latitude'
        lon_name = 'lon' if 'lon' in data.dims else 'longitude'
        weights = np.cos(np.deg2rad(data[lat_name]))
        avg = data.weighted(weights).mean(dim=[lat_name, lon_name]).values
        tp_list.append(avg)
    return np.array(tp_list)

# ===== Get trends for each region =====
nw_tp = read_weighted_tp(paths['NW'], 'prev')
ne_tp = read_weighted_tp(paths['NE'], 'prev')
sw_tp = read_weighted_tp(paths['SW'], 'prev')
se_tp = read_weighted_tp(paths['SE'], 'prev')

slope_nw, *_ = linregress(years, nw_tp)
slope_ne, *_ = linregress(years, ne_tp)
slope_sw, *_ = linregress(years, sw_tp)
slope_se, *_ = linregress(years, se_tp)

# ===== Relative contribution calculation =====
external_total = slope_nw + slope_ne + slope_sw + slope_se
contributions = {
    'NW': slope_nw / external_total,
    'NE': slope_ne / external_total,
    'SW': slope_sw / external_total,
    'SE': slope_se / external_total
}

# ===== Draw bar chart =====
labels = list(contributions.keys())
sizes = [v * 100 for v in contributions.values()]  # Convert to percentage
colors = ['#CA0E12', '#25377F', '#2AA7DE', '#F6BD21']

x = np.arange(len(labels))  # X-axis positions
bar_width = 0.4

ax_f = fig.add_subplot(3,2,6)
bars = ax_f.bar(x, sizes, width=bar_width, color=colors, edgecolor='black')

# Add percentage text labels
for bar, size in zip(bars, sizes):
    ax_f.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
              f'{size:.1f}%', ha='center', va='bottom', fontsize=20)

# Graphic beautification
ax_f.set_xticks(x)
ax_f.set_xticklabels(labels, fontsize=16)
ax_f.tick_params(axis='y', labelsize=16)
ax_f.set_ylabel('Percentage (%)', fontsize=18)
ax_f.set_title('Relative contribution to trend', fontsize=20)
ax_f.set_ylim(0, max(sizes) * 1.15)
ax_f.text(0.00, 1.02, "(f)", transform=ax_f.transAxes, fontsize=20, ha='left', va='bottom', fontweight='bold')

plt.subplots_adjust(hspace=0.25, wspace=0.25)
plt.savefig('',dpi =800, bbox_inches="tight")