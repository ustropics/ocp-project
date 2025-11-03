#!/usr/bin/env python3
# --------------------------------------------------------------
#  Depth-averaged Salinity Transport |S × V|
#  → Three-segment rectangles per band WITH GAPS:
#       • 32–30°N: FILLED (alpha=0.8), black border
#       • 30–20°N: BORDER ONLY (no fill), colored edge
#       • 20–18°N: FILLED (alpha=0.8), black border
#     → Small gaps (0.08°) between segments to show all borders
# --------------------------------------------------------------

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime
import os
from netCDF4 import Dataset
from matplotlib.patches import Rectangle

# ------------------------------------------------------------------
# 1. File paths
# ------------------------------------------------------------------
cmap = 'cubehelix'

file_path = f'output/data/25deg/3d/gomb4_daily_2024_{day}_3z.nc'
save_dir = 'output/images'
os.makedirs(save_dir, exist_ok=True)
save_path = os.path.join(save_dir, f'salinity_transport_2024_{day}_bands.png')

# ------------------------------------------------------------------
# 2. FOUR LONGITUDE BANDS
# ------------------------------------------------------------------
bands = [
    {"west": -96.0, "east": -93.0, "label": "96°W–93°W"},
    {"west": -93.0, "east": -90.0, "label": "93°W–90°W"},
    {"west": -90.0, "east": -87.0, "label": "90°W–87°W"},
    {"west": -87.0, "east": -84.0, "label": "87°W–84°W"},
]

# ------------------------------------------------------------------
# 3. HARD-CODED BAND COLOURS (matches plt_salinx_bands.py)
# ------------------------------------------------------------------
BAND_COLORS = [
    (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),   # blue
    (1.0,                 0.4980392156862745, 0.054901960784313725),  # orange
    (0.17254901960784313, 0.6274509803921569, 0.17254901960784313),  # green
    (0.8392156862745098,  0.15294117647058825, 0.1568627450980392),  # red
]

# ------------------------------------------------------------------
# 4. LATITUDE SEGMENTS + GAP
# ------------------------------------------------------------------
gap = 0.08  # degrees (approximately 9 km) — adjust if needed

raw_segments = [
    {"north": 32.0, "south": 30.0, "fill": True,  "alpha": 0.8},
    {"north": 30.0, "south": 20.0, "fill": False, "alpha": None},
    {"north": 20.0, "south": 18.0, "fill": True,  "alpha": 0.8},
]

# Apply gap: shrink each segment inward by gap/2 on both ends
segments = []
for i, seg in enumerate(raw_segments):
    north = seg["north"]
    south = seg["south"]

    # Apply gap only between segments (not at outer edges)
    if i > 0:
        south += gap / 2
    if i < len(raw_segments) - 1:
        north -= gap / 2

    segments.append({
        "north": north,
        "south": south,
        "fill": seg["fill"],
        "alpha": seg["alpha"]
    })

# ------------------------------------------------------------------
# 5. Get VALID_RANGE
# ------------------------------------------------------------------
nc = Dataset(file_path, 'r')
valid_s = nc.variables['salinity'].valid_range
valid_u = nc.variables['u'].valid_range
valid_v = nc.variables['v'].valid_range
nc.close()

# ------------------------------------------------------------------
# 6. Load data
# ------------------------------------------------------------------
ds = xr.open_dataset(file_path)
S = ds['salinity'].squeeze()
u = ds['u'].squeeze()
v = ds['v'].squeeze()
lat = ds['Latitude'].values
lon = ds['Longitude'].values

# ------------------------------------------------------------------
# 7. Timestamp
# ------------------------------------------------------------------
date_val = ds['Date'].values.item()
date_str = f"{int(date_val):08d}"
date_obj = datetime.strptime(date_str, "%Y%m%d")
month_name = date_obj.strftime("%B")
day_num = date_obj.day


def ordinal(n):
    return "%d%s" % (n, "tsnrhtdd"[(n // 10 % 10 != 1) * (n % 10 < 4) * n % 10::4])


timestamp = f"{month_name} {ordinal(day_num)}, {date_obj.year}"

# ------------------------------------------------------------------
# 8. MASK & compute transport
# ------------------------------------------------------------------
S = S.where((S >= valid_s[0]) & (S <= valid_s[1]), np.nan)
u = u.where((u >= valid_u[0]) & (u <= valid_u[1]), np.nan)
v = v.where((v >= valid_v[0]) & (v <= valid_v[1]), np.nan)

Su = S * u
Sv = S * v
Su_mean = Su.mean(dim='Depth', skipna=True)
Sv_mean = Sv.mean(dim='Depth', skipna=True)
transport_mag = np.sqrt(Su_mean**2 + Sv_mean**2)
transport_mag = transport_mag.where(np.isfinite(transport_mag), np.nan)

# ------------------------------------------------------------------
# 9. Stats
# ------------------------------------------------------------------
mag_valid = transport_mag.stack(z=('Latitude', 'Longitude')).dropna('z')
mag_min = float(mag_valid.min()) if len(mag_valid) > 0 else np.nan
mag_max = float(mag_valid.max()) if len(mag_valid) > 0 else np.nan
mag_mean = float(mag_valid.mean()) if len(mag_valid) > 0 else np.nan

# ------------------------------------------------------------------
# 10. Plot setup
# ------------------------------------------------------------------
fig = plt.figure(figsize=(11, 8))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([lon.min(), lon.max(), lat.min(),
              lat.max()], crs=ccrs.PlateCarree())

# Contourf
levels = np.linspace(0, 40, 41)
cf = ax.contourf(lon, lat, transport_mag, levels=levels,
                 cmap=cmap, extend='max', transform=ccrs.PlateCarree())

# Contours
cs = ax.contour(lon, lat, transport_mag, levels=np.linspace(0, 30, 4),
                colors='white', linewidths=0.5, transform=ccrs.PlateCarree())
cs2 = ax.contour(lon, lat, transport_mag, levels=np.linspace(30, 40, 2),
                 colors='black', linewidths=0.5, transform=ccrs.PlateCarree())
ax.clabel(cs, inline=True, fontsize=9, fmt='%.2f', colors='white')
ax.clabel(cs2, inline=True, fontsize=9, fmt='%.2f', colors='black')

# Colorbar
cbar = plt.colorbar(cf, ax=ax, shrink=0.7, pad=0.06)
cbar.set_label('Salinity Transport Magnitude |S × V| (psu·m/s)', fontsize=12)

# Map features
ax.coastlines(resolution='10m', linewidth=1, color='black', zorder=3)
ax.add_feature(cfeature.LAND, facecolor="#9c6b00", zorder=2)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, zorder=3, color='black')
ax.add_feature(cfeature.LAKES, facecolor='#a0c8f0',
               edgecolor='black', linewidth=0.5, zorder=3)
ax.add_feature(cfeature.STATES, linewidth=0.5, zorder=3,
               facecolor='none', edgecolor='black')

gl = ax.gridlines(draw_labels=True, linestyle='--', color='gray', alpha=0.5)
gl.top_labels = gl.right_labels = False

ax.set_title(
    f'Depth-Averaged Salinity Transport ({timestamp})', fontsize=14, pad=15)
# ax.text(0.01, 0.01, timestamp, transform=ax.transAxes, fontsize=11, color='white',
#         bbox=dict(facecolor='black', alpha=0.7, edgecolor='none', pad=3))

# ------------------------------------------------------------------
# 11. DRAW THREE-SEGMENT RECTANGLES WITH GAPS
# ------------------------------------------------------------------
for i, b in enumerate(bands):
    color = BAND_COLORS[i]
    west, east = b["west"], b["east"]
    width = east - west

    for j, seg in enumerate(segments):
        north, south = seg["north"], seg["south"]
        height = north - south

        if seg["fill"]:
            # Top & bottom: filled, strong alpha, black border
            facecolor = color
            edgecolor = 'black'
            alpha = seg["alpha"]
            linewidth = 1.2
        else:
            # Middle: border only, colored edge, no fill
            facecolor = color
            edgecolor = color
            alpha = .25
            linewidth = 1  # slightly thicker for visibility

        rect = Rectangle(
            (west, south), width, height,
            linewidth=linewidth,
            edgecolor=edgecolor,
            facecolor=facecolor,
            alpha=alpha,
            transform=ccrs.PlateCarree(),
            zorder=4
        )
        ax.add_patch(rect)

        # Label only in middle (transparent) segment
        if j == 1:
            label_x = (west + east) / 2.0
            label_y = (north + south) / 1.59
            ax.text(label_x, label_y, b["label"],
                    ha='center', va='bottom', fontsize=8.5, color='black',
                    transform=ccrs.PlateCarree(),
                    bbox=dict(facecolor='white', alpha=0.95,
                              edgecolor='none', pad=1.5),
                    zorder=5)

# ------------------------------------------------------------------
# 12. Save
# ------------------------------------------------------------------
plt.tight_layout()
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.close(fig)
print(f"Plot saved: {save_path}")
