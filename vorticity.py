import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import os
import ssl
import geopandas as gpd
from metpy.calc import vorticity, lat_lon_grid_deltas
from metpy.units import units
import pandas as pd

ssl._create_default_https_context = ssl._create_unverified_context

# File path to the NetCDF file
file_path = 'era5/era5_file.nc'

# Load the NetCDF file using xarray
ds = xr.open_dataset(file_path)

# Define the domain (latitudes: 4N to 22N, longitudes: 114E to 130E)
lat_min, lat_max = 4, 22
lon_min, lon_max = 114, 130

# Filter the dataset by the domain
ds = ds.sel(latitude=slice(lat_max, lat_min), longitude=slice(lon_min, lon_max))

# Define output directory
output_dir = 'wind_vort_output'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Load the land shapefile
shapefile = 'shapefiles/phprov.shp'
land = gpd.read_file(shapefile)

# Function to sanitize filenames
def sanitize_filename(filename):
    return filename.replace(":", "-").replace(".", "-")

# Loop through all valid times and pressure levels
for time_idx in range(len(ds.valid_time)):
    if ds.valid_time[time_idx].dt.hour.values != 0:
        continue

    for level_idx in range(len(ds.pressure_level)):
        # Extract wind components and lat/lon
        u = ds.u.isel(valid_time=time_idx, pressure_level=level_idx).metpy.sel(
            latitude=slice(lat_max, lat_min), longitude=slice(lon_min, lon_max)).squeeze().values * units('m/s')
        v = ds.v.isel(valid_time=time_idx, pressure_level=level_idx).metpy.sel(
            latitude=slice(lat_max, lat_min), longitude=slice(lon_min, lon_max)).squeeze().values * units('m/s')
        
        lat = ds.latitude.values
        lon = ds.longitude.values

        # Compute dx and dy in meters using MetPy helper function
        dx, dy = lat_lon_grid_deltas(lon, lat)

        # Compute relative vorticity
        vort = vorticity(u, v, dx=dx, dy=dy)

        # Create output directory for this pressure level
        level_name = f"{int(ds.pressure_level[level_idx].values)} mb"
        level_dir = os.path.join(output_dir, level_name)
        if not os.path.exists(level_dir):
            os.makedirs(level_dir)

        # Plotting
        plt.figure(figsize=(7, 7))
        c = plt.pcolormesh(lon, lat, vort.magnitude, cmap='bwr', shading='auto', vmin=-5e-5, vmax=5e-5)
        cbar = plt.colorbar(c, label="Relative Vorticity (s⁻¹)")
        cbar.set_ticks(np.arange(-5e-5, 5.1e-5, 1e-5))

        valid_time = pd.to_datetime(ds.valid_time[time_idx].values)
        valid_time_str = valid_time.strftime('%Y-%m-%d-%H')
        plt.title(f"{level_name} Relative Vorticity - {valid_time_str} UTC")

        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        land.boundary.plot(ax=plt.gca(), color='black', linewidth=1)

        # Save the plot
        plot_filename = f"vorticity_{int(ds.pressure_level[level_idx].values)}_{valid_time_str}.png"
        plot_filepath = os.path.join(level_dir, plot_filename)

        plt.savefig(plot_filepath)
        plt.close()

        print(f"Finished: {plot_filename}")

ds.close()
print(f"Relative vorticity charts saved in {output_dir}")
