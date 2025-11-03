#!/usr/bin/env python3
"""
plt_salin_multi_band_cross.py
Latitudinal cross-sections (22°–28°N) of depth-averaged |S × V|
for 4 specified longitude bands.
Shows AVERAGE (not peak) per band in stats box.
→ Uses HARD-CODED colours to match plt_salin.py
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import os
from netCDF4 import Dataset

# ------------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------------
file_path = f'output/data/25deg/3d/gomb4_daily_2024_{day}_3z.nc'
save_dir = 'output/images'
os.makedirs(save_dir, exist_ok=True)
save_path = os.path.join(
    save_dir, f'salinity_transport_multi_band_2024_{day}.png')

# === FOUR BANDS ===
bands = [
    {"west": -96.0, "east": -93.0, "label": "96°W–93°W"},
    {"west": -93.0, "east": -90.0, "label": "93°W–90°W"},
    {"west": -90.0, "east": -87.0, "label": "90°W–87°W"},
    {"west": -87.0, "east": -84.0, "label": "87°W–84°W"},
]

# Latitude limit
lat_min, lat_max = 22.0, 28.0

# ------------------------------------------------------------------
# 2. HARD-CODED BAND COLOURS (identical to plt_salin.py)
# ------------------------------------------------------------------
BAND_COLORS = [
    (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),   # blue
    (1.0,                 0.4980392156862745, 0.054901960784313725),  # orange
    (0.17254901960784313, 0.6274509803921569, 0.17254901960784313),  # green
    (0.8392156862745098,  0.15294117647058825, 0.1568627450980392),  # red
]

# ------------------------------------------------------------------
# 3. Load valid_range
# ------------------------------------------------------------------
nc = Dataset(file_path, 'r')
valid_s = nc.variables['salinity'].valid_range
valid_u = nc.variables['u'].valid_range
valid_v = nc.variables['v'].valid_range
nc.close()

# ------------------------------------------------------------------
# 4. Load data + latitude slice
# ------------------------------------------------------------------
ds = xr.open_dataset(file_path)
ds = ds.sel(Latitude=slice(lat_min, lat_max))

S = ds['salinity'].squeeze()
u = ds['u'].squeeze()
v = ds['v'].squeeze()
lat = ds['Latitude'].values
lon = ds['Longitude'].values

# ------------------------------------------------------------------
# 5. Timestamp
# ------------------------------------------------------------------
date_val = ds['Date'].values.item()
date_str = f"{int(date_val):08d}"
date_obj = datetime.strptime(date_str, "%Y%m%d")
month_name = date_obj.strftime("%B")


def ordinal(n):
    return "%d%s" % (n, "tsnrhtdd"[(n//10 % 10 != 1)*(n % 10 < 4)*n % 10::4])


timestamp = f"{month_name} {ordinal(date_obj.day)}, {date_obj.year}"

# ------------------------------------------------------------------
# 6. Mask & compute transport
# ------------------------------------------------------------------
S = S.where((S >= valid_s[0]) & (S <= valid_s[1]), np.nan)
u = u.where((u >= valid_u[0]) & (u <= valid_u[1]), np.nan)
v = v.where((v >= valid_v[0]) & (v <= valid_v[1]), np.nan)

Su = S * u
Sv = S * v
transport_mag = np.sqrt(Su**2 + Sv**2)
transport_avg = transport_mag.mean(dim='Depth', skipna=True)

# ------------------------------------------------------------------
# 7. Zonal mean for each band
# ------------------------------------------------------------------
profiles = []
for b in bands:
    mask = (lon >= b["west"]) & (lon <= b["east"])
    band_data = transport_avg.where(mask)
    zonal = band_data.mean(dim='Longitude', skipna=True)
    profiles.append(zonal)

# ------------------------------------------------------------------
# 8. Plot
# ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(11, 7))

for i, (prof, b) in enumerate(zip(profiles, bands)):
    ax.plot(lat, prof, label=b["label"], color=BAND_COLORS[i], linewidth=2.8)
    ax.fill_between(lat, prof, 0, color=BAND_COLORS[i], alpha=0.2)

ax.set_xlabel('Latitude (°N)', fontsize=13)
ax.set_ylabel('Depth-Averaged |S × V| (psu·m/s)', fontsize=13)
ax.set_title(f'Four-Band Salinity Transport Cross-Sections (22°–28°N)\n{timestamp}',
             fontsize=14, pad=15)

ax.grid(True, linestyle='--', alpha=0.6)
ax.legend(fontsize=11, title="Longitude Bands",
          title_fontsize=12, loc='upper left')
ax.set_xlim(lat_min, lat_max)

# ------------------------------------------------------------------
# 9. Stats: AVERAGE per band
# ------------------------------------------------------------------
stats_lines = []
for i, b in enumerate(bands):
    avg_val = profiles[i].mean(skipna=True).item()
    stats_lines.append(f"{b['label']}: {avg_val:.3f}")

ax.text(0.98, 0.98, f"Avg per band:\n" + "\n".join(stats_lines),
        transform=ax.transAxes, fontsize=10,
        verticalalignment='top', horizontalalignment='right',
        bbox=dict(facecolor='white', alpha=0.95, boxstyle="round,pad=0.5", edgecolor='gray'))

# ------------------------------------------------------------------
# 10. Save
# ------------------------------------------------------------------
plt.tight_layout()
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.close(fig)
print(f"Multi-band plot saved: {save_path}")
