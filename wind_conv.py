import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import os
import ssl
import geopandas as gpd
from metpy.calc import divergence
from metpy.units import units
import pandas as pd  # Import pandas for date handling
from matplotlib.ticker import MaxNLocator  # Import for granular color bar ticks

ssl._create_default_https_context = ssl._create_unverified_context

# File path to the NetCDF file
file_path = 'era5/8c426675800d082614e8c529e1d8b0c1.nc'

# Load the NetCDF file using xarray
ds = xr.open_dataset(file_path)

# Define the domain (latitudes: 4N to 22N, longitudes: 114E to 130E)
lat_min, lat_max = 4, 22
lon_min, lon_max = 114, 130

# Filter the dataset by the domain
ds = ds.sel(latitude=slice(lat_max, lat_min), longitude=slice(lon_min, lon_max))

# Define output directory
output_dir = 'wind_conv_output'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Load the land shapefile (shapefile location changed to 'shapefiles/phprov.shp')
shapefile = 'shapefiles/phprov.shp'
land = gpd.read_file(shapefile)

# Function to compute convergence (using MetPy's divergence function)
def compute_convergence(u, v, lat, lon):
    """
    Compute the wind convergence (or divergence) using the u and v wind components.
    Convergence is the negative of the divergence.
    """
    # Manually compute gradients using central difference
    # Compute the gradient in the longitude (x) direction
    du_dx = (np.roll(u, -1, axis=-1) - u) / (lon[1] - lon[0])
    
    # Compute the gradient in the latitude (y) direction
    dv_dy = (np.roll(v, -1, axis=-2) - v) / (lat[1] - lat[0])

    # Calculate divergence using the gradients
    div = du_dx + dv_dy

    #returns negative of divergence (which is convergence)
    return -div

# Function to sanitize filenames
def sanitize_filename(filename):
    return filename.replace(":", "-").replace(".", "-")

# Loop through all valid times and pressure levels, then plot convergence
for time_idx in range(len(ds.valid_time)):
    # Only process 0000 UTC times
    if ds.valid_time[time_idx].dt.hour.values != 0:
        continue
    
    for level_idx in range(len(ds.pressure_level)):
        
        # Extract wind components for the current valid time and pressure level
        u_data = ds.u.isel(valid_time=time_idx, pressure_level=level_idx).metpy.sel(latitude=slice(lat_max, lat_min), longitude=slice(lon_min, lon_max))
        v_data = ds.v.isel(valid_time=time_idx, pressure_level=level_idx).metpy.sel(latitude=slice(lat_max, lat_min), longitude=slice(lon_min, lon_max))

        # Extract the latitude and longitude arrays
        lat = ds.latitude.values
        lon = ds.longitude.values
        
        # Ensure the data are 2D arrays (latitude, longitude)
        u_data = u_data.squeeze()  # Remove any singleton dimensions
        v_data = v_data.squeeze()  # Remove any singleton dimensions

        # Calculate convergence (divergence)
        div = compute_convergence(u_data, v_data, lat, lon)
        
        # Create the output folder for the current level if not already existing
        level_name = f"{int(ds.pressure_level[level_idx].values)} mb"  # Convert to integer to remove decimals
        level_dir = os.path.join(output_dir, level_name)
        if not os.path.exists(level_dir):
            os.makedirs(level_dir)
        
        # Create a plot for convergence (shaded plot)
        plt.figure(figsize=(7, 7))
        # Use pcolormesh for shaded areas with finer levels
        c = plt.pcolormesh(lon, lat, div, cmap='coolwarm', shading='auto', vmin=-0.00025, vmax=0.00025)
        
        # Create the color bar and set its ticks
        cbar = plt.colorbar(c, label="Convergence (s^-1)")
        cbar.set_ticks(np.arange(-0.00025, 0.00025, 0.00005))  # Set finer tick intervals

        # Modify the title format as requested
        valid_time = pd.to_datetime(ds.valid_time[time_idx].values)
        valid_time_str = valid_time.strftime('%Y-%m-%d-%H')  # Format as 'yyyy-mm-dd-HH'
        plt.title(f"{int(ds.pressure_level[level_idx].values)} mb Wind Convergence - {valid_time_str} UTC")
        
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        
        # Overlay the land shapefile
        land.boundary.plot(ax=plt.gca(), color='black', linewidth=1)
        
        # Sanitize the filename to avoid invalid characters
        sanitized_time = sanitize_filename(str(ds.valid_time[time_idx].values))
        
        # Create the filename with the desired format: convergence_(level)_yyyy-mm-dd-tt.png
        plot_filename = f"convergence_{int(ds.pressure_level[level_idx].values)}_{valid_time_str}.png"
        plot_filepath = os.path.join(level_dir, plot_filename)
        
        # Save the plot as a .png file in the appropriate directory
        plt.savefig(plot_filepath)
        plt.close()

        # Print that the file is done
        print(f"Finished: {plot_filename}")

# Close the dataset
ds.close()

print(f"Convergence charts saved in {output_dir}")
