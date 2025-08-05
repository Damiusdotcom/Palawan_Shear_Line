# Script created by DJEV – Updated for domain control, zorder, and tidy tick formatting
import os
import ssl
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import geopandas as gpd
import cartopy.crs as ccrs
import cartopy.feature as cfeature

ssl._create_default_https_context = ssl._create_unverified_context

# === Configuration ===
csv_file = "rainfall_data.csv"
shapefile_path = "shapefiles/phprov.shp"
output_dir = "rainfall_actual"
os.makedirs(output_dir, exist_ok=True)

# === Domain / Map extent ===
domain = [113, 128, 5, 21]  # [lon_min, lon_max, lat_min, lat_max]

# === Load data ===
df = pd.read_csv(csv_file)
df.columns = df.columns.str.strip()
date_columns = [col for col in df.columns if col.isdigit() and len(col) == 8]

df["lat"] = pd.to_numeric(df["lat"], errors="coerce")
df["lon"] = pd.to_numeric(df["lon"], errors="coerce")
shapefile = gpd.read_file(shapefile_path)

# === Define precipitation colormap ===
levels = [0, 1, 10, 25, 50, 100, 200, 300, 500, 700]
colors = ["#ffffff", "#bab8b8", "#00c5ff", "#6bfb90", "#ffff00", "#ffaa00", "#ff0000", "#ff73df", "#8400a8"]
cmap = mcolors.ListedColormap(colors).with_extremes(over="#000000")
norm = mcolors.BoundaryNorm(levels, cmap.N)

# === Plotting loop ===
for date in date_columns:
    df[date] = df[date].replace("T", 0)
    df[date] = pd.to_numeric(df[date], errors="coerce")

    fig, ax = plt.subplots(figsize=(8, 6), subplot_kw={"projection": ccrs.PlateCarree()})
    ax.set_extent(domain, crs=ccrs.PlateCarree())

    # Map features
    ax.add_feature(cfeature.BORDERS, linestyle=":")
    ax.add_feature(cfeature.COASTLINE)
    shapefile.boundary.plot(ax=ax, edgecolor="black", linewidth=1, transform=ccrs.PlateCarree())

    # Scatter rainfall with zorder
    sc = ax.scatter(
        df["lon"], df["lat"],
        c=df[date], cmap=cmap, norm=norm,
        edgecolor="black", s=50, zorder=3,
        transform=ccrs.PlateCarree()
    )

    # Colorbar
    cbar = plt.colorbar(sc, ax=ax, orientation="vertical", shrink=0.7, pad=0.02, extend="max")
    cbar.set_label("Precipitation (mm)")
    cbar.set_ticks(levels)
    cbar.set_ticklabels([str(l) for l in levels])

    # Title with formatted date
    formatted_date = pd.to_datetime(date, format="%Y%m%d").strftime("%Y-%m-%d")
    ax.set_title(f"Synoptic Station 24-hour Accumulated Rainfall", fontsize=14)

    # Lat/lon ticks (multiples of 5)
    lon_ticks = np.arange(np.ceil(domain[0] / 5) * 5, domain[1] + 1, 5)
    lat_ticks = np.arange(np.ceil(domain[2] / 5) * 5, domain[3] + 1, 5)

    ax.set_xticks(lon_ticks, crs=ccrs.PlateCarree())
    ax.set_yticks(lat_ticks, crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, pos: f"{v:.0f}°E"))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, pos: f"{v:.0f}°N"))
    ax.tick_params(axis="both", labelsize=10)
    ax.set_xlabel("")
    ax.set_ylabel("")

    # Save
    plt.savefig(f"{output_dir}/{date}.png", dpi=300, bbox_inches="tight")
    plt.close()

print(f"Plots saved in '{output_dir}'")
