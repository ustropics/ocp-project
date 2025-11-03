#!/usr/bin/env python3
"""
plt_salinx_vert_combined.py
4 individual + 1 combined 2x2 plot
Inset map at [0.345, 0.1, 0.28, 0.28] in individuals
Inset maps in each subplot of combined figure
Band highlighted in BAND_COLORS
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import os
from netCDF4 import Dataset
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# ------------------------------------------------------------------
# 1. CONFIG
# ------------------------------------------------------------------
cmap = 'cubehelix'
file_path = f'output/data/25deg/3d/gomb4_daily_2024_{day}_3z.nc'
save_dir = 'output/images'
os.makedirs(save_dir, exist_ok=True)

# FOUR BANDS + COLORS
BAND_COLORS = [
    (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),   # blue
    (1.0,                 0.4980392156862745, 0.054901960784313725),  # orange
    (0.17254901960784313, 0.6274509803921569, 0.17254901960784313),  # green
    (0.8392156862745098,  0.15294117647058825, 0.1568627450980392),  # red
]

bands = [
    {"west": -96.0, "east": -93.0, "label": "96°W–93°W",
        "slug": "96W-93W", "color": BAND_COLORS[0]},
    {"west": -93.0, "east": -90.0, "label": "93°W–90°W",
        "slug": "93W-90W", "color": BAND_COLORS[1]},
    {"west": -90.0, "east": -87.0, "label": "90°W–87°W",
        "slug": "90W-87W", "color": BAND_COLORS[2]},
    {"west": -87.0, "east": -84.0, "label": "87°W–84°W",
        "slug": "87W-84W", "color": BAND_COLORS[3]},
]

lat_min, lat_max = 22.0, 28.0
depth_min, depth_max = 0.0, 1000.0

# ------------------------------------------------------------------
# 2. LOAD VALID RANGE
# ------------------------------------------------------------------
nc = Dataset(file_path, 'r')
valid_s = nc.variables['salinity'].valid_range
valid_u = nc.variables['u'].valid_range
valid_v = nc.variables['v'].valid_range
nc.close()

# ------------------------------------------------------------------
# 3. LOAD DATA
# ------------------------------------------------------------------
ds_full = xr.open_dataset(file_path)
ds = ds_full.sel(Latitude=slice(lat_min, lat_max),
                 Depth=slice(depth_min, depth_max))

S = ds['salinity'].squeeze()
u = ds['u'].squeeze()
v = ds['v'].squeeze()
lat = ds['Latitude'].values
lon = ds['Longitude'].values
depth = ds['Depth'].values

# Full domain for inset
lon_full = ds_full['Longitude'].values
lat_full = ds_full['Latitude'].values

# ------------------------------------------------------------------
# 4. TIMESTAMP
# ------------------------------------------------------------------
date_val = int(ds_full['Date'].values.item())
date_obj = datetime.strptime(f"{date_val:08d}", "%Y%m%d")
month_name = date_obj.strftime("%B")


def ordinal(n): return "%d%s" % (
    n, "tsnrhtdd"[(n//10 % 10 != 1)*(n % 10 < 4)*n % 10::4])


timestamp = f"{month_name} {ordinal(date_obj.day)}, {date_obj.year}"

# ------------------------------------------------------------------
# 5. MASK
# ------------------------------------------------------------------
S = S.where((S >= valid_s[0]) & (S <= valid_s[1]), np.nan)
u = u.where((u >= valid_u[0]) & (u <= valid_u[1]), np.nan)
v = v.where((v >= valid_v[0]) & (v <= valid_v[1]), np.nan)

# ------------------------------------------------------------------
# 6. TRANSPORT
# ------------------------------------------------------------------
transport_mag = np.sqrt((S * u)**2 + (S * v)**2)

# ------------------------------------------------------------------
# 7. GLOBAL COLOR SCALE
# ------------------------------------------------------------------
all_zonal = []
for b in bands:
    mask = (lon >= b["west"]) & (lon <= b["east"])
    zonal = transport_mag.where(mask).mean(dim='Longitude', skipna=True)
    all_zonal.append(zonal.values)
flat = np.concatenate([z.ravel() for z in all_zonal if np.any(~np.isnan(z))])
vmax = np.nanpercentile(flat, 98)
levels = np.linspace(0, max(40, vmax), 41)

# ------------------------------------------------------------------
# 8. PRE-COMPUTE ZONAL DATA
# ------------------------------------------------------------------
zonal_data = []
means = []
for b in bands:
    mask = (lon >= b["west"]) & (lon <= b["east"])
    zonal = transport_mag.where(mask).mean(dim='Longitude', skipna=True)
    zonal_data.append(zonal.values)
    means.append(zonal.mean(skipna=True).item())

# ------------------------------------------------------------------
# 9. INDIVIDUAL PLOTS WITH INSET MAP (at [0.345, 0.1, 0.28, 0.28])
# ------------------------------------------------------------------
print("\nGenerating 4 individual plots with inset maps...\n")
for i, b in enumerate(bands):
    west, east = b["west"], b["east"]
    label, slug = b["label"], b["slug"]
    band_color = b["color"]
    save_path = os.path.join(
        save_dir, f'salinity_transport_vertical_{slug}_2024_{day}.png')

    Z = zonal_data[i]
    X, Y = np.meshgrid(lat, depth)

    # Main plot
    fig = plt.figure(figsize=(11, 7))
    ax_main = fig.add_axes([0.1, 0.1, 0.65, 0.8])

    cf = ax_main.contourf(X, Y, Z, levels=levels, cmap=cmap, extend='max')
    cs = ax_main.contour(
        X, Y, Z, levels=levels[::5], colors='white', linewidths=0.6)
    ax_main.clabel(cs, inline=True, fontsize=9, fmt='%.1f', colors='white')

    cbar = fig.colorbar(cf, ax=ax_main, shrink=0.8, pad=0.05)
    cbar.set_label('|S × V| (psu·m/s)', fontsize=12)

    ax_main.set_xlabel('Latitude (°N)', fontsize=13)
    ax_main.set_ylabel('Depth (m)', fontsize=13)
    ax_main.set_title(
        f'Vertical Cross‑Section\n{label} | {timestamp}', fontsize=14, pad=15)
    ax_main.set_ylim(depth_max, depth_min)
    ax_main.grid(True, linestyle='--', alpha=0.5)

    ax_main.text(0.02, 0.02, f"0–1000 m avg:\n{means[i]:.3f} psu·m/s",
                 transform=ax_main.transAxes, fontsize=10, ha='left', va='bottom',
                 bbox=dict(facecolor='white', alpha=0.9, boxstyle="round,pad=0.4", edgecolor='gray'))

    # === INSET MAP at [0.345, 0.1, 0.28, 0.28] ===
    ax_inset = fig.add_axes([0.424, 0.1, 0.2, 0.2],
                            projection=ccrs.PlateCarree())
    ax_inset.set_extent([lon_full.min(), lon_full.max(
    ), lat_full.min(), lat_full.max()], crs=ccrs.PlateCarree())

    ax_inset.coastlines(resolution='10m', linewidth=0.8)
    ax_inset.add_feature(cfeature.LAND, facecolor='#d9d9d9')
    ax_inset.add_feature(cfeature.BORDERS, linewidth=0.5)

    # Highlight current band
    rect = Rectangle((west, lat_min), east - west, lat_max - lat_min,
                     linewidth=2.5, edgecolor=band_color, facecolor=band_color, alpha=0.6,
                     transform=ccrs.PlateCarree())
    ax_inset.add_patch(rect)

    # Latitude slice
    ax_inset.hlines([lat_min, lat_max], lon_full.min(), lon_full.max(),
                    color='black', linewidth=1.5, linestyle='--', transform=ccrs.PlateCarree())

    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"SAVED: {save_path}")

# ------------------------------------------------------------------
# 10. COMBINED 2x2 PLOT (with inset maps in every panel)
# ------------------------------------------------------------------
print("\nGenerating combined 2x2 plot with inset maps in each panel...\n")

fig, axes = plt.subplots(2, 2, figsize=(15, 11), sharex=True, sharey=True)
axes = axes.flatten()

# ---- inset tuning -------------------------------------------------------
inset_rel_w = 0.30          # width  of inset as fraction of subplot
inset_rel_h = 0.30          # height of inset as fraction of subplot
inset_pad_x = 0.08          # base right-shift inside subplot
inset_pad_y = 0.08          # base down-shift inside subplot

# Pixel-to-figure conversions (300 dpi)
PIXELS_DOWN_GLOBAL = 46 / (fig.dpi * fig.get_figheight())   # ~0.015
PIXELS_UP_TOP = 23 / (fig.dpi * fig.get_figheight())   # net 20 px down for top
# +10 px right (left panels)
PIXELS_RIGHT_LEFT = 6 / (fig.dpi * fig.get_figwidth())
# +20 px right (right panels)
PIXELS_RIGHT_RIGHT = 21 / (fig.dpi * fig.get_figwidth())

for i, (b, Z, mean_val) in enumerate(zip(bands, zonal_data, means)):
    ax = axes[i]
    X, Y = np.meshgrid(lat, depth)

    # ---- main contour plot ---------------------------------------------
    cf = ax.contourf(X, Y, Z, levels=levels, cmap=cmap, extend='max')
    cs = ax.contour(
        X, Y, Z, levels=levels[::5], colors='white', linewidths=0.6)
    ax.clabel(cs, inline=True, fontsize=8, fmt='%.1f', colors='white')

    ax.set_title(b["label"], fontsize=13, pad=8,
                 fontweight='bold', color=b["color"])
    ax.set_ylim(depth_max, depth_min)
    ax.grid(True, linestyle='--', alpha=0.4)

    ax.text(0.98, 0.98, f"{mean_val:.3f}", transform=ax.transAxes,
            fontsize=9, color='white', fontweight='bold',
            ha='right', va='top',
            bbox=dict(facecolor='black', alpha=0.7, pad=2))

    # ---- inset map (lower-right of *this* subplot) --------------------
    bbox = ax.get_position()

    # base position (lower-right)
    inset_left = bbox.x0 + (1 - inset_rel_w - inset_pad_x) * bbox.width
    inset_bottom = bbox.y0 + inset_pad_y * bbox.height

    # **Global 50 px down**
    inset_bottom -= PIXELS_DOWN_GLOBAL

    # **Top row: move up total 30 px → net 20 px down**
    if i < 2:
        inset_bottom += PIXELS_UP_TOP

    # **Horizontal shift: left vs right panels**
    if i % 2 == 0:  # left column (i = 0, 2)
        inset_left += PIXELS_RIGHT_LEFT
    else:           # right column (i = 1, 3)
        inset_left += 0.04 * bbox.width + PIXELS_RIGHT_RIGHT

    # final size
    inset_width = inset_rel_w * bbox.width
    inset_height = inset_rel_h * bbox.height

    ax_inset = fig.add_axes(
        [inset_left, inset_bottom, inset_width, inset_height],
        projection=ccrs.PlateCarree()
    )

    ax_inset.set_extent([lon_full.min(), lon_full.max(),
                         lat_full.min(), lat_full.max()], crs=ccrs.PlateCarree())

    ax_inset.coastlines(resolution='10m', linewidth=0.6)
    ax_inset.add_feature(cfeature.LAND, facecolor='#d9d9d9', alpha=0.8)
    ax_inset.add_feature(cfeature.BORDERS, linewidth=0.4)

    # highlight current band
    rect = Rectangle((b["west"], lat_min), b["east"] - b["west"],
                     lat_max - lat_min,
                     linewidth=2, edgecolor=b["color"],
                     facecolor=b["color"], alpha=0.6,
                     transform=ccrs.PlateCarree())
    ax_inset.add_patch(rect)

    # latitude slice
    ax_inset.hlines([lat_min, lat_max], lon_full.min(), lon_full.max(),
                    color='black', linewidth=1.2, linestyle='--',
                    transform=ccrs.PlateCarree())

    ax_inset.set_xticks([])
    ax_inset.set_yticks([])

# ---- shared labels & colorbar -------------------------------------------
fig.supxlabel('Latitude (°N)', fontsize=14, y=0.04)
fig.supylabel('Depth (m)', fontsize=14, x=0.04)

cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
cbar = fig.colorbar(cf, cax=cbar_ax)
cbar.set_label('|S × V| (psu·m/s)', fontsize=12)

fig.suptitle(f'Vertical Cross-Sections of Salinity Transport (0–1000 m)\n{timestamp}',
             fontsize=16, y=0.95)

combined_path = os.path.join(save_dir,
                             f'salinity_transport_vertical_4panel_2024_{day}.png')
plt.subplots_adjust(left=0.08, right=0.90, top=0.88,
                    bottom=0.10, wspace=0.20, hspace=0.30)
plt.savefig(combined_path, dpi=300, bbox_inches='tight')
plt.close(fig)
print(f"SAVED COMBINED: {combined_path}")

print("\nAll 5 plots (4 individual + 1 combined with insets) generated successfully!")
