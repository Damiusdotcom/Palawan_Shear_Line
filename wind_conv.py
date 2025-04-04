import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import os

# File path to the NetCDF file
file_path = 'era5/8c426675800d082614e8c529e1d8b0c1.nc'

# Load the NetCDF file using xarray
ds = xr.open_dataset(file_path)

# Define the domain (latitudes: 0 to 27N, longitudes: 110E to 155E)
lat_min, lat_max = 0, 27
lon_min, lon_max = 110, 155

# Filter the dataset by the domain
ds = ds.sel(latitude=slice(lat_max, lat_min), longitude=slice(lon_min, lon_max))

# Define output directory
output_dir = 'wind_conv_output'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Function to compute convergence (simple divergence calculation)
def compute_convergence(u, v, dx, dy):
    """
    Compute the wind convergence (or divergence) using the u and v wind components.
    Convergence is the negative of the divergence.
    """
    # Ensure u and v are numpy arrays (no metadata)
    u = np.asarray(u)  # Use np.asarray to convert to ndarray if it's not already
    v = np.asarray(v)  # Convert to ndarray

    # Print the data type and shape of u and v for debugging
    print(f"Shape of u: {u.shape}, dtype: {u.dtype}")
    print(f"Shape of v: {v.shape}, dtype: {v.dtype}")

    # Compute derivatives (central difference) using numpy.roll
    du_dx = (np.roll(u, shift=1, axis=-1) - u) / dx  # Shift in longitude (x direction)
    dv_dy = (np.roll(v, shift=1, axis=-2) - v) / dy  # Shift in latitude (y direction)
    
    # Convergence = - (du/dx + dv/dy)
    convergence = -(du_dx + dv_dy)

    # Ensure that the output is 2D (latitude, longitude)
    convergence = np.squeeze(convergence)  # Remove any singleton dimensions
    return convergence

# Get the latitude and longitude resolution (assuming regular grid)
dx = ds.longitude[1] - ds.longitude[0]
dy = ds.latitude[1] - ds.latitude[0]

# Loop through all valid times and pressure levels, then plot convergence
for time_idx in range(len(ds.valid_time)):
    for level_idx in range(len(ds.pressure_level)):
        
        # Extract wind components for the current valid time and pressure level
        u_data = ds.u.isel(valid_time=time_idx, pressure_level=level_idx).values
        v_data = ds.v.isel(valid_time=time_idx, pressure_level=level_idx).values
        
        # Print the dimensions of u_data and v_data to debug
        print(f"Extracted u_data shape: {u_data.shape}")
        print(f"Extracted v_data shape: {v_data.shape}")
        
        # Ensure the data are 2D arrays (latitude, longitude)
        u_data = np.squeeze(u_data)  # Remove any singleton dimensions
        v_data = np.squeeze(v_data)  # Remove any singleton dimensions
        
        # Calculate convergence
        convergence = compute_convergence(u_data, v_data, dx, dy)
        
        # Create the output folder for the current level if not already existing
        level_name = f"{ds.pressure_level[level_idx].values} mb"
        level_dir = os.path.join(output_dir, level_name)
        if not os.path.exists(level_dir):
            os.makedirs(level_dir)
        
        # Create a plot for convergence
        plt.figure(figsize=(12, 8))
        plt.contourf(ds.longitude, ds.latitude, convergence, cmap='coolwarm', levels=np.linspace(-0.0005, 0.0005, 20))
        plt.colorbar(label="Convergence (s^-1)")
        plt.title(f"Wind Convergence at {str(ds.valid_time[time_idx].values)} for {level_name}")
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        
        # Save the plot as a .png file in the appropriate directory
        plot_filename = f"convergence_{str(ds.valid_time[time_idx].values)}.png"
        plot_filepath = os.path.join(level_dir, plot_filename)
        plt.savefig(plot_filepath)
        plt.close()

# Close the dataset
ds.close()

print(f"Convergence charts saved in {output_dir}")
