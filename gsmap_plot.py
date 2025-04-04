import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
import matplotlib.colors as mcolors
import shapefile as shp
from cartopy.io.shapereader import Reader

# Define input and output directories
data_dir = "gsmap_data"
output_dir = "gsmap_plots"
os.makedirs(output_dir, exist_ok=True)

# Define domain for map limits (not slicing the dataset)
lat_min, lat_max = 4, 22   # Latitude range
lon_min, lon_max = 114, 130  # Longitude range

# Define precipitation thresholds (mm) and corresponding colors
levels = [0, 1, 10, 25, 50, 100, 200, 300, 500, 700]
colors = ["#ffffff", "#bab8b8", "#00c5ff", "#6bfb90", "#ffff00", "#ffaa00", "#ff0000", "#ff73df", "#8400a8"]

# Create a custom colormap and norm
cmap = mcolors.ListedColormap(colors).with_extremes(over="#000000")
norm = mcolors.BoundaryNorm(levels, cmap.N)

# Loop through each subfolder (each day)
for day_folder in sorted(os.listdir(data_dir)):
    day_path = os.path.join(data_dir, day_folder)
    if not os.path.isdir(day_path):
        continue

    # Get all NetCDF files for the day
    nc_files = sorted([os.path.join(day_path, f) for f in os.listdir(day_path) if f.endswith(".nc")])

    if not nc_files:
        print(f"No NetCDF files found for {day_folder}")
        continue

    accumulated_rainfall = None
    lon, lat = None, None

    for nc_file in nc_files:
        try:
            dataset = xr.open_dataset(nc_file)
        except ValueError:
            dataset = xr.open_dataset(nc_file, decode_times=False)

        # Extract precipitation variable (hourly precipitation rate)
        if "hourlyPrecipRateGC" in dataset:
            rain = dataset["hourlyPrecipRateGC"].squeeze()
        else:
            print(f"Skipping {nc_file}, 'hourlyPrecipRateGC' not found.")
            dataset.close()
            continue

        # Ensure dataset is not empty
        if rain.size == 0:
            print(f"Skipping {nc_file}, data is empty.")
            dataset.close()
            continue

        # Store latitude and longitude only from the first valid dataset
        if lon is None and lat is None:
            lon = dataset.Longitude.values
            lat = dataset.Latitude.values

        # Ensure rain matches expected dimensions
        if accumulated_rainfall is None:
            accumulated_rainfall = np.zeros_like(rain)

        accumulated_rainfall += rain.fillna(0)  # Replace NaNs with 0
        dataset.close()

    # Ensure accumulated_rainfall is not empty
    if accumulated_rainfall is None or np.all(accumulated_rainfall == 0):
        print(f"No valid rainfall data for {day_folder}")
        continue

    # Ensure accumulated_rainfall has the correct shape
    accumulated_rainfall = accumulated_rainfall.squeeze()

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={"projection": ccrs.PlateCarree()})
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())  # Set map limits

    # Add map features
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray')
    
    # Overlay shapefile
    shapefile_path = "shapefiles/phprov.shp"
    ax.add_geometries(Reader(shapefile_path).geometries(), ccrs.PlateCarree(), edgecolor='black', facecolor='none')

    # Plot rainfall
    cf = ax.pcolormesh(lon, lat, accumulated_rainfall, cmap=cmap, norm=norm, transform=ccrs.PlateCarree(), shading='auto')

    # Add a custom color bar
    cbar = plt.colorbar(cf, ax=ax, orientation="vertical", shrink=0.7, pad=0.02, extend="max")
    cbar.set_label("Precipitation (mm)")
    cbar.set_ticks(levels)
    cbar.set_ticklabels([f"{int(l)}" for l in levels])

    # Add latitude and longitude bars with 5-degree intervals
    ax.set_xticks(np.arange(lon_min, lon_max + 1, 5), crs=ccrs.PlateCarree())
    ax.set_yticks(np.arange(lat_min, lat_max + 1, 5), crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: f"{val}\u00b0E"))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: f"{val}\u00b0N"))

    # Add title
    ax.set_title(f"GSMAP 24-hour Accumulated Rainfall {day_folder[:4]}-{day_folder[4:6]}-{day_folder[6:]}")

    # Save plot
    output_filename = os.path.join(output_dir, f"{day_folder}.png")
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.close(fig)

    print(f"Finished outputting: {output_filename}")

print("All 24-hour accumulated rainfall images successfully generated and saved.")
