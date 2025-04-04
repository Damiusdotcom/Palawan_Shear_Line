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

# Ensure output directory exists
output_dir = "gsmap_point"
os.makedirs(output_dir, exist_ok=True)

# Load CSV file with coordinates (assumed columns: 'lat', 'lon')
csv_file = "coordinates.csv"
coords_df = pd.read_csv(csv_file)

# Define region boundaries (adjust as needed)
lat_min, lat_max = 4, 22
lon_min, lon_max = 114, 130

# Load the shapefile
shapefile = gpd.read_file("shapefiles/phprov.shp")

# Define precipitation thresholds and colors
levels = [0, 1, 10, 25, 50, 100, 200, 300, 500, 700]
colors = ["#ffffff", "#bab8b8", "#00c5ff", "#6bfb90", "#ffff00", "#ffaa00", "#ff0000", "#ff73df", "#8400a8"]
cmap = mcolors.ListedColormap(colors).with_extremes(over="#000000")
norm = mcolors.BoundaryNorm(levels, cmap.N)

# Loop through all NetCDF files in the 'gsmap_data' folder
for folder in sorted(glob.glob("gsmap_data/*")):
    folder_name = os.path.basename(folder)
    print(f"Processing: {folder_name}")

    # Accumulate 24-hour rainfall from all available NetCDF files
    gsmap_total = None
    nc_files = glob.glob(os.path.join(folder, "*.nc"))
    if not nc_files:
        print(f"Skipping {folder_name}, no GSMaP NetCDF files found.")
        continue
    
    for nc_file in nc_files:
        ds = xr.open_dataset(nc_file, decode_times=False)
        
        # Extract rainfall variable (adjust variable name if needed)
        rainfall = ds["hourlyPrecipRate"].values.squeeze()  # Correct variable name
        lats = ds["Latitude"].values  # Use 'Latitude' instead of 'lat'
        lons = ds["Longitude"].values  # Use 'Longitude' instead of 'lon'
        
        if gsmap_total is None:
            gsmap_total = rainfall
        else:
            gsmap_total += rainfall
        
        ds.close()
    
    if gsmap_total is None:
        print(f"Skipping {folder_name}, no valid GSMaP data found.")
        continue
    
    # Extract date from folder name (yyyymmdd)
    date_str = folder_name

    # Extract precipitation data at each coordinate
    precip_values = []
    for _, row in coords_df.iterrows():
        lat, lon = row["lat"], row["lon"]
        precip_value = gsmap_total[np.abs(lats - lat).argmin(), np.abs(lons - lon).argmin()]
        precip_values.append(precip_value)

    # Add precipitation data to DataFrame
    coords_df["precip"] = precip_values

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6), subplot_kw={"projection": ccrs.PlateCarree()})
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    # Add map features
    ax.add_feature(cfeature.BORDERS, linestyle=":")
    ax.add_feature(cfeature.COASTLINE)

    # Plot shapefile
    shapefile.boundary.plot(ax=ax, edgecolor="black", linewidth=1, transform=ccrs.PlateCarree())

    # Scatter plot with precipitation data
    sc = ax.scatter(
        coords_df["lon"], coords_df["lat"], 
        c=coords_df["precip"], cmap=cmap, norm=norm,
        edgecolor="black", s=50, transform=ccrs.PlateCarree()
    )

    # Define tick values for longitude and latitude at 5-degree increments
    lon_ticks = np.arange(lon_min, lon_max + 1, 5)
    lat_ticks = np.arange(lat_min, lat_max + 1, 5)

    # Set tick positions
    ax.set_xticks(lon_ticks, crs=ccrs.PlateCarree())
    ax.set_yticks(lat_ticks, crs=ccrs.PlateCarree())

    # Format tick labels to display degrees
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, pos: f"{v:.1f}°E"))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, pos: f"{v:.1f}°N"))

    # Adjust tick label size
    ax.tick_params(axis="both", labelsize=10)

    # Remove labels
    ax.set_xlabel("")
    ax.set_ylabel("")

    # Add color bar
    cbar = plt.colorbar(sc, ax=ax, orientation="vertical", shrink=0.7, pad=0.02, extend="max")
    cbar.set_label("Precipitation (mm)")
    cbar.set_ticks(levels)
    cbar.set_ticklabels([str(l) for l in levels])

    # Set title with extracted date in 'yyyy-mm-dd' format
    start_date = datetime.strptime(date_str, "%Y%m%d")
    ax.set_title(f"GSMap 24-hour Accumulated Rainfall ({start_date.strftime('%Y-%m-%d')})", fontsize=14)

    # Save plot
    plt.savefig(f"gsmap_point/{start_date.strftime('%Y%m%d')}.png", dpi=300, bbox_inches='tight')
    plt.close()
    
print("Processing complete!")
