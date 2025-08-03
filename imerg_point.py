import os
import glob
import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import geopandas as gpd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime

# === CONFIGURATION ===
input_dir = "imerg_finalrun"
output_dir = "imerg_point_output"
csv_file = "coordinates.csv"
shapefile_path = "shapefiles/phprov.shp"
os.makedirs(output_dir, exist_ok=True)

# Region bounds
lat_min, lat_max = 4, 22
lon_min, lon_max = 114, 130

# Color thresholds and map
levels = [0, 1, 10, 25, 50, 100, 200, 300, 500, 700]
colors = ["#ffffff", "#bab8b8", "#00c5ff", "#6bfb90", "#ffff00", "#ffaa00", "#ff0000", "#ff73df", "#8400a8"]
cmap = mcolors.ListedColormap(colors).with_extremes(over="#000000")
norm = mcolors.BoundaryNorm(levels, cmap.N)

# Load coordinate CSV
coords_df = pd.read_csv(csv_file)

# Load shapefile
shapefile = gpd.read_file(shapefile_path)


def plot_imerg_points(nc_path):
    print(f"Processing: {nc_path}")
    dataset = xr.open_dataset(nc_path)

    # Extract date from NetCDF
    date_str = pd.to_datetime(str(dataset.time[0].values)).strftime("%Y-%m-%d")
    date_file = pd.to_datetime(str(dataset.time[0].values)).strftime("%Y%m%d")

    # Extract precipitation at each coordinate
    precip_values = []
    for _, row in coords_df.iterrows():
        lat, lon = row["lat"], row["lon"]
        precip_value = dataset["precipitation"].sel(lat=lat, lon=lon, method="nearest").values
        precip_values.append(precip_value)

    # Add precipitation data to DataFrame
    coords_df["precip"] = precip_values

    # Create plot
    fig, ax = plt.subplots(figsize=(8, 6), subplot_kw={"projection": ccrs.PlateCarree()})
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    ax.add_feature(cfeature.BORDERS, linestyle=":")
    ax.add_feature(cfeature.COASTLINE)
    shapefile.boundary.plot(ax=ax, edgecolor="black", linewidth=1, transform=ccrs.PlateCarree())

    sc = ax.scatter(
        coords_df["lon"], coords_df["lat"],
        c=coords_df["precip"], cmap=cmap, norm=norm,
        edgecolor="black", s=50, transform=ccrs.PlateCarree()
    )

    # Set ticks and labels
    lon_ticks = np.arange(lon_min, lon_max + 1, 5)
    lat_ticks = np.arange(lat_min, lat_max + 1, 5)
    ax.set_xticks(lon_ticks, crs=ccrs.PlateCarree())
    ax.set_yticks(lat_ticks, crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, pos: f"{v:.1f}°E"))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, pos: f"{v:.1f}°N"))
    ax.tick_params(axis="both", labelsize=10)
    ax.set_xlabel("")
    ax.set_ylabel("")

    cbar = plt.colorbar(sc, ax=ax, orientation="vertical", shrink=0.7, pad=0.02, extend="max")
    cbar.set_label("Precipitation (mm)")
    cbar.set_ticks(levels)
    cbar.set_ticklabels([str(l) for l in levels])

    ax.set_title(f"IMERG 24-hour Accumulated Rainfall ({date_str})", fontsize=14)

    plt.savefig(os.path.join(output_dir, f"{date_file}.png"), dpi=300, bbox_inches="tight")
    plt.close()


# === MAIN LOOP ===
for nc_file in sorted(glob.glob(os.path.join(input_dir, "*.nc4"))):
    plot_imerg_points(nc_file)

print("Processing complete!")
