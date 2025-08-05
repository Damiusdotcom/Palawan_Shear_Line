import os
import glob
import ssl
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
from shapely.geometry import box
from datetime import datetime

ssl._create_default_https_context = ssl._create_unverified_context

# === CONFIGURATION ===
input_dir = "imerg_finalrun"
output_dir = "imerg_plot_output"
shapefile_path = "shapefiles/phprov.shp"
os.makedirs(output_dir, exist_ok=True)

# Map region boundaries
lat_min, lat_max = 5, 21
lon_min, lon_max = 113, 128

# Color thresholds and colormap
levels = [0, 1, 10, 25, 50, 100, 200, 300, 500, 700]
colors = ["#e0e0e0", "#bab8b8", "#00c5ff", "#6bfb90", "#ffff00", "#ffaa00", "#ff0000", "#ff73df", "#8400a8"]
cmap = mcolors.ListedColormap(colors).with_extremes(over="#000000")
norm = mcolors.BoundaryNorm(levels, cmap.N)

# Load and clip shapefile
shapefile = gpd.read_file(shapefile_path)
bbox = box(lon_min, lat_min, lon_max, lat_max)
shapefile = gpd.clip(shapefile, bbox)


def plot_precipitation(nc_path):
    dataset = xr.open_dataset(nc_path)
    precip = dataset["precipitation"]
    regional_precip = precip.sel(lat=slice(lat_min, lat_max), lon=slice(lon_min, lon_max))

    # Sum over time if more than one step
    precip_slice = (
        regional_precip.sum(dim="time") if "time" in regional_precip.dims and regional_precip.sizes["time"] > 1
        else regional_precip.isel(time=0)
    )

    # Extract date for title and filename
    date = pd.to_datetime(str(dataset.time[0].values)).strftime("%Y-%m-%d")
    date_file = pd.to_datetime(str(dataset.time[0].values)).strftime("%Y%m%d")

    # Create figure and axes
    fig, ax = plt.subplots(figsize=(6, 6), subplot_kw={"projection": ccrs.PlateCarree()})
    im = precip_slice.plot.imshow(
        ax=ax, cmap=cmap, norm=norm, x="lon", y="lat", transform=ccrs.PlateCarree(), add_colorbar=False
    )

    shapefile.boundary.plot(ax=ax, edgecolor="black", linewidth=1, transform=ccrs.PlateCarree())
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=":")

    ax.set_title(f"IMERG 24-hour Accumulated Rainfall", fontsize=14)

    # Assuming lon_min, lon_max, lat_min, lat_max are already defined
    lon_ticks = np.arange(
        ((lon_min + 4) // 5) * 5,  # Round up to nearest multiple of 5
        ((lon_max) // 5) * 5 + 1,  # Round down to nearest multiple of 5
        5
    )
    lat_ticks = np.arange(
        ((lat_min + 4) // 5) * 5, 
        ((lat_max) // 5) * 5 + 1, 
        5
    )

    ax.set_xticks(lon_ticks, crs=ccrs.PlateCarree())
    ax.set_yticks(lat_ticks, crs=ccrs.PlateCarree())

    ax.xaxis.set_major_formatter(plt.FuncFormatter(
        lambda v, pos: f"{abs(v):.0f}°{'E' if v >= 0 else 'W'}"
    ))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(
        lambda v, pos: f"{abs(v):.0f}°{'N' if v >= 0 else 'S'}"
    ))

    ax.tick_params(axis="both", labelsize=10)
    ax.set_xlabel("")
    ax.set_ylabel("")

    cbar = plt.colorbar(im, ax=ax, orientation="vertical", shrink=0.7, pad=0.02, extend="max")
    cbar.set_label("Precipitation (mm)")
    cbar.set_ticks(levels)
    cbar.set_ticklabels([f"{int(l)}" for l in levels])

    plt.savefig(os.path.join(output_dir, f"{date_file}.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Processed: {nc_path}")


# === MAIN LOOP ===
for nc_file in sorted(glob.glob(os.path.join(input_dir, "*.nc4"))):
    plot_precipitation(nc_file)
