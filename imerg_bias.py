import os
import glob
import ssl
import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime

ssl._create_default_https_context = ssl._create_unverified_context

# === CONFIGURATION ===
input_dir = "imerg"
output_dir = "imerg_bias_output"
station_file = "rainfall_data.csv"
shapefile_path = "shapefiles/phprov.shp"
os.makedirs(output_dir, exist_ok=True)

# Load station data
station_df = pd.read_csv(station_file)
station_df.columns = station_df.columns.str.strip()
station_df["lat"] = pd.to_numeric(station_df["lat"], errors="coerce")
station_df["lon"] = pd.to_numeric(station_df["lon"], errors="coerce")
rainfall_dates = [col for col in station_df.columns if col.isdigit() and len(col) == 8]

# Load shapefile
shapefile = gpd.read_file(shapefile_path)

# Color scale and bias levels
levels = [-100, -50, -25, -10, -5, -1, 0, 1, 5, 10, 25, 50, 100]
colors = ["#a80000", "#ff0000", "#ff7373", "#ffb8b8", "#ffe6e6", "#ffffff",
          "#e6f7ff", "#b8e6ff", "#73c2ff", "#7373ff", "#0000a8", "#000066"]
cmap = mcolors.ListedColormap(colors)
norm = mcolors.BoundaryNorm(levels, cmap.N)


def plot_imerg_station_bias(imerg_path):
    dataset = xr.open_dataset(imerg_path)
    date_str = pd.to_datetime(str(dataset.time[0].values)).strftime("%Y-%m-%d")
    date_col = date_str.replace("-", "")

    if date_col not in rainfall_dates:
        print(f"Skipping {date_str}, no matching station data.")
        return

    # Convert and clean station values
    station_df[date_col] = station_df[date_col].replace("T", 0).astype(float)

    # Match IMERG to station coordinates
    imerg_values = [
        dataset["precipitation"].sel(lat=row["lat"], lon=row["lon"], method="nearest").values
        for _, row in station_df.iterrows()
    ]

    # Calculate bias
    station_df["imerg"] = imerg_values
    station_df["bias"] = station_df["imerg"] - station_df[date_col]

    # Create plot
    fig, ax = plt.subplots(figsize=(8, 6), subplot_kw={"projection": ccrs.PlateCarree()})
    ax.set_extent([
        station_df["lon"].min() - 1, station_df["lon"].max() + 1,
        station_df["lat"].min() - 1, station_df["lat"].max() + 1
    ], crs=ccrs.PlateCarree())

    ax.add_feature(cfeature.BORDERS, linestyle=":")
    ax.add_feature(cfeature.COASTLINE)
    shapefile.boundary.plot(ax=ax, edgecolor="black", linewidth=1, transform=ccrs.PlateCarree())

    # Gridlines
    gl = ax.gridlines(draw_labels=True, linestyle="--", linewidth=0.5)
    gl.top_labels = False
    gl.right_labels = False

    # Scatter plot
    sc = ax.scatter(
        station_df["lon"], station_df["lat"],
        c=station_df["bias"], cmap=cmap, norm=norm,
        edgecolor="black", s=50, transform=ccrs.PlateCarree()
    )

    # Colorbar
    cbar = plt.colorbar(sc, ax=ax, orientation="vertical", shrink=0.7, pad=0.02, extend="both")
    cbar.set_label("Bias (IMERG - Station) mm")
    cbar.set_ticks(levels)
    cbar.set_ticklabels([str(l) for l in levels])

    ax.set_title(f"(IMERG - Station) Bias ({date_str})", fontsize=14)
    plt.savefig(os.path.join(output_dir, f"{date_col}.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Processed: {date_str}")


# === MAIN LOOP ===
for file in sorted(glob.glob(os.path.join(input_dir, "*.nc4"))):
    plot_imerg_station_bias(file)

print("Bias plots saved in 'imerg_bias_output'")
