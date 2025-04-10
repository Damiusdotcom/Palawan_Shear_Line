import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import os
import ssl
import geopandas as gpd
from metpy.units import units
from metpy.calc import lat_lon_grid_deltas, divergence
import pandas as pd

ssl._create_default_https_context = ssl._create_unverified_context

# File path to the NetCDF file
file_path = 'era5/era5_file.nc'
ds = xr.open_dataset(file_path)

# Domain
lat_min, lat_max = 4, 22
lon_min, lon_max = 114, 130
ds = ds.sel(latitude=slice(lat_max, lat_min), longitude=slice(lon_min, lon_max))

# Output directory
output_dir = 'moisture_flux_conv_output'
os.makedirs(output_dir, exist_ok=True)

# Load shapefile
shapefile = 'shapefiles/phprov.shp'
land = gpd.read_file(shapefile)

def sanitize_filename(filename):
    return filename.replace(":", "-").replace(".", "-")

# Loop through all 0000 UTC valid times
for time_idx in range(len(ds.valid_time)):
    if ds.valid_time[time_idx].dt.hour.values != 0:
        continue

    for level_idx in range(len(ds.pressure_level)):
        level_val = int(ds.pressure_level[level_idx].values)
        level_name = f"{level_val} mb"
        level_dir = os.path.join(output_dir, level_name)
        os.makedirs(level_dir, exist_ok=True)

        # Extract data
        q = ds.q.isel(valid_time=time_idx, pressure_level=level_idx).metpy.sel(
            latitude=slice(lat_max, lat_min), longitude=slice(lon_min, lon_max)).squeeze().values * units('kg/kg')
        u = ds.u.isel(valid_time=time_idx, pressure_level=level_idx).metpy.sel(
            latitude=slice(lat_max, lat_min), longitude=slice(lon_min, lon_max)).squeeze().values * units('m/s')
        v = ds.v.isel(valid_time=time_idx, pressure_level=level_idx).metpy.sel(
            latitude=slice(lat_max, lat_min), longitude=slice(lon_min, lon_max)).squeeze().values * units('m/s')

        lat = ds.latitude.values
        lon = ds.longitude.values
        dx, dy = lat_lon_grid_deltas(lon, lat)

        # Moisture flux components
        qu = q * u
        qv = q * v

        # Compute divergence of moisture flux (negative for convergence)
        mfc = -divergence(qu, qv, dx=dx, dy=dy)

        # Plotting
        plt.figure(figsize=(7, 7))
        c = plt.pcolormesh(lon, lat, mfc.magnitude, cmap='BrBG', shading='auto', vmin=-5e-6, vmax=5e-6)
        cbar = plt.colorbar(c, label='Moisture Flux Convergence (s⁻¹)')
        cbar.set_ticks(np.arange(-5e-6, 5.1e-6, 1e-6))

        valid_time = pd.to_datetime(ds.valid_time[time_idx].values)
        valid_time_str = valid_time.strftime('%Y-%m-%d-%H')
        plt.title(f"{level_name} Moisture Flux Convergence - {valid_time_str} UTC")

        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
        land.boundary.plot(ax=plt.gca(), color='black', linewidth=1)

        # Save plot
        plot_filename = f"mfc_{level_val}_{valid_time_str}.png"
        plot_filepath = os.path.join(level_dir, plot_filename)
        plt.savefig(plot_filepath)
        plt.close()

        print(f"Finished: {plot_filename}")

ds.close()
print(f"Moisture flux convergence charts saved in {output_dir}")
