import os
import glob
import xarray as xr
import ssl
import matplotlib.pyplot as plt
import geopandas as gpd
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from shapely.geometry import box
import numpy as np
from datetime import datetime, timedelta

ssl._create_default_https_context = ssl._create_unverified_context

# Define the region
lat_min, lat_max = 4, 22  # Adjust as per your region
lon_min, lon_max = 114, 130

# Load the shapefile (replace with correct path)
shapefile = gpd.read_file("shapefiles/phprov.shp")

# Clip shapefile using bounding box
bbox = box(lon_min, lat_min, lon_max, lat_max)
shapefile = gpd.clip(shapefile, bbox)

# Define precipitation thresholds (mm) and corresponding colors
levels = [0, 1, 10, 25, 50, 100, 200, 300, 500, 700]
colors = ["#ffffff", "#bab8b8", "#00c5ff", "#6bfb90", "#ffff00", "#ffaa00", "#ff0000", "#ff73df", "#8400a8"]

# Create a custom colormap and norm
cmap = mcolors.ListedColormap(colors).with_extremes(over="#000000")
norm = mcolors.BoundaryNorm(levels, cmap.N)

# Loop through all NetCDF files in the imerg directory
for file in sorted(glob.glob("imerg/*.nc4")):
    dataset = xr.open_dataset(file)
    
    # Extract precipitation data
    precip = dataset["precipitation"]
    regional_precip = precip.sel(lat=slice(lat_min, lat_max), lon=slice(lon_min, lon_max))
    precip_slice = regional_precip.isel(time=0)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(6, 6), subplot_kw={"projection": ccrs.PlateCarree()})
    
    # Plot precipitation with thresholds
    im = precip_slice.plot.imshow(ax=ax, cmap=cmap, norm=norm, x="lon", y="lat", transform=ccrs.PlateCarree(), add_colorbar=False)
    
    # Plot shapefile boundaries 
    shapefile.boundary.plot(ax=ax, edgecolor="black", linewidth=1, transform=ccrs.PlateCarree())
    
    # Restrict map extent
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    
    # Add map features
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=":")
    
    # Extract date from the dataset
    date_str = str(dataset.time[0].values)[:10]  # Extract YYYY-MM-DD from datetime format

    # Set title
    ax.set_title(f"IMERG 24-hour Accumulated Rainfall {date_str}", fontsize=14)

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
    
    # Add a custom color bar
    cbar = plt.colorbar(im, ax=ax, orientation="vertical", shrink=0.7, pad=0.02, extend="max")
    cbar.set_label("Precipitation (mm)")
    cbar.set_ticks(levels)
    cbar.set_ticklabels([f"{int(l)}" for l in levels])
    
    # Save plot
    start_date = np.datetime64(dataset.time[0].values)  # Extract first time value
    start_date = datetime.utcfromtimestamp(start_date.astype('O') / 1e9)  # Convert to Python datetime
    start_date_str = start_date.strftime("%B_%d_%Y")  # Format as "Month_Day_Year"

    plt.savefig(f"imerg_plots/{start_date.strftime('%Y%m%d')}.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Processed: {file}")
