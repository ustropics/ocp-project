#!/usr/bin/env python3
"""
plt_salin_faceted.py
Latitudinal cross-sections (22°–28°N) of depth-averaged |S × V|
for 4 specified longitude bands.
Each band is shown in its own subplot (small multiples).
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
    save_dir, f'salinity_transport_faceted_2024_{day}.png')

# === FOUR BANDS ===
bands = [
    {"west": -96.0, "east": -93.0, "label": "96°W–93°W"},
    {"west": -93.0, "east": -90.0, "label": "93°W–90°W"},
    {"west": -90.0, "east": -87.0, "label": "90°W–87°W"},
    {"west": -87.0, "east": -84.0, "label": "87°W–84°W"},
]

lat_min, lat_max = 22.0, 28.0

BAND_COLORS = [
    (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
    (1.0, 0.4980392156862745, 0.054901960784313725),
    (0.17254901960784313, 0.6274509803921569, 0.17254901960784313),
    (0.8392156862745098, 0.15294117647058825, 0.1568627450980392),
]

# ------------------------------------------------------------------
# 2. Load valid_range
# ------------------------------------------------------------------
nc = Dataset(file_path, 'r')
valid_s = nc.variables['salinity'].valid_range
valid_u = nc.variables['u'].valid_range
valid_v = nc.variables['v'].valid_range
nc.close()

# ------------------------------------------------------------------
# 3. Load data
# ------------------------------------------------------------------
ds = xr.open_dataset(file_path)
ds = ds.sel(Latitude=slice(lat_min, lat_max))
S = ds['salinity'].squeeze()
u = ds['u'].squeeze()
v = ds['v'].squeeze()
lat = ds['Latitude'].values
lon = ds['Longitude'].values

# ------------------------------------------------------------------
# 4. Timestamp
# ------------------------------------------------------------------
date_val = ds['Date'].values.item()
date_str = f"{int(date_val):08d}"
date_obj = datetime.strptime(date_str, "%Y%m%d")
month_name = date_obj.strftime("%B")


def ordinal(n): return "%d%s" % (
    n, "tsnrhtdd"[(n//10 % 10 != 1)*(n % 10 < 4)*n % 10::4])


timestamp = f"{month_name} {ordinal(date_obj.day)}, {date_obj.year}"

# ------------------------------------------------------------------
# 5. Mask & compute transport
# ------------------------------------------------------------------
S = S.where((S >= valid_s[0]) & (S <= valid_s[1]), np.nan)
u = u.where((u >= valid_u[0]) & (u <= valid_u[1]), np.nan)
v = v.where((v >= valid_v[0]) & (v <= valid_v[1]), np.nan)

Su = S * u
Sv = S * v
transport_mag = np.sqrt(Su**2 + Sv**2)
transport_avg = transport_mag.mean(dim='Depth', skipna=True)

# ------------------------------------------------------------------
# 6. Zonal mean for each band
# ------------------------------------------------------------------
profiles = []
for b in bands:
    mask = (lon >= b["west"]) & (lon <= b["east"])
    zonal = transport_avg.where(mask).mean(dim='Longitude', skipna=True)
    profiles.append(zonal)

# ------------------------------------------------------------------
# 7. Faceted plot
# ------------------------------------------------------------------
fig, axes = plt.subplots(4, 1, figsize=(10, 10), sharex=True, sharey=True)

for i, (prof, b) in enumerate(zip(profiles, bands)):
    ax = axes[i]
    ax.plot(lat, prof, color=BAND_COLORS[i], linewidth=2.5)
    ax.fill_between(lat, prof, 0, color=BAND_COLORS[i], alpha=0.25)
    ax.set_title(b["label"], fontsize=11, loc='left')
    ax.grid(True, linestyle='--', alpha=0.5)

axes[-1].set_xlabel('Latitude (°N)')
axes[1].set_ylabel('Depth-Averaged |S × V| (psu·m/s)')

fig.suptitle(f'Salinity Transport by Longitude Band\n{timestamp}', fontsize=14)
plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig(save_path, dpi=300)
plt.close(fig)
print(f"Faceted plot saved: {save_path}")
