import os
import glob
import pandas as pd
import numpy as np
import ssl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import geopandas as gpd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
from datetime import datetime

ssl._create_default_https_context = ssl._create_unverified_context

# Ensure output directory exists
output_dir = "bias_imerg-station"
os.makedirs(output_dir, exist_ok=True)

# Load station rainfall data
station_file = "rainfall_data.csv"
station_df = pd.read_csv(station_file)
station_df.columns = station_df.columns.str.strip()
station_df["lat"] = pd.to_numeric(station_df["lat"], errors="coerce")
station_df["lon"] = pd.to_numeric(station_df["lon"], errors="coerce")

# Identify rainfall columns dynamically
rainfall_dates = [col for col in station_df.columns if col.isdigit() and len(col) == 8]

# Load the shapefile
shapefile = gpd.read_file("shapefiles/phprov.shp")

# Define bias color scale
levels = [-100, -50, -25, -10, -5, -1, 0, 1, 5, 10, 25, 50, 100]
colors = ["#a80000", "#ff0000", "#ff7373", "#ffb8b8", "#ffe6e6", "#ffffff", "#e6f7ff", "#b8e6ff", "#73c2ff", "#7373ff", "#0000a8", "#000066"]
cmap = mcolors.ListedColormap(colors)
norm = mcolors.BoundaryNorm(levels, cmap.N)

# Loop through each NetCDF file
for imerg_file in sorted(glob.glob("imerg/*.nc4")):
    dataset = xr.open_dataset(imerg_file)
    date_str = str(dataset.time[0].values)[:10]  # Extract YYYY-MM-DD
    date_col = date_str.replace("-", "")  # Convert to YYYYMMDD format
    
    if date_col not in rainfall_dates:
        print(f"Skipping {date_str}, no matching station data.")
        continue
    
    # Get station rainfall data as reference
    station_df[date_col] = station_df[date_col].replace("T", 0).astype(float)
    
    # Get IMERG rainfall at station locations
    imerg_values = []
    for _, row in station_df.iterrows():
        lat, lon = row["lat"], row["lon"]
        precip_value = dataset["precipitation"].sel(lat=lat, lon=lon, method="nearest").values
        imerg_values.append(precip_value)
    
    station_df["imerg"] = imerg_values
    
    # Calculate bias (IMERG - Station)
    station_df["bias"] = station_df["imerg"] - station_df[date_col]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6), subplot_kw={"projection": ccrs.PlateCarree()})
    
    # Set map extent
    ax.set_extent([station_df["lon"].min() - 1, station_df["lon"].max() + 1, 
                   station_df["lat"].min() - 1, station_df["lat"].max() + 1], crs=ccrs.PlateCarree())
    
    # Add map features
    ax.add_feature(cfeature.BORDERS, linestyle=":")
    ax.add_feature(cfeature.COASTLINE)
    shapefile.boundary.plot(ax=ax, edgecolor="black", linewidth=1, transform=ccrs.PlateCarree())
    
    # Add latitude and longitude grid lines
    gl = ax.gridlines(draw_labels=True, linestyle="--", linewidth=0.5)
    gl.top_labels = False
    gl.right_labels = False
    
    # Scatter plot with bias data
    sc = ax.scatter(station_df["lon"], station_df["lat"], c=station_df["bias"], cmap=cmap, norm=norm,
                     edgecolor="black", s=50, transform=ccrs.PlateCarree())
    
    # Add color bar
    cbar = plt.colorbar(sc, ax=ax, orientation="vertical", shrink=0.7, pad=0.02, extend="both")
    cbar.set_label("Bias (IMERG - Station) mm")
    cbar.set_ticks(levels)
    cbar.set_ticklabels([str(l) for l in levels])
    
    # Set title
    ax.set_title(f"(IMERG-Station) Bias ({date_str})", fontsize=14)
    
    # Save plot
    plt.savefig(f"{output_dir}/{date_col}.png", dpi=300, bbox_inches="tight")
    plt.close()
    
    print(f"Processed: {date_str}")

print("Bias plots saved in 'bias_imerg-station'")