import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
import matplotlib.colors as mcolors

# Define file path
nc_file = "era5/era5_file.nc"

# Open dataset
dataset = xr.open_dataset(nc_file)

# Define domain (Modify these values as needed)
lat_min, lat_max = 0, 27   # Latitude range
lon_min, lon_max = 110, 145  # Longitude range

# Identify available pressure levels in dataset
pressure_levels = dataset.pressure_level.values.tolist()

# Output directory
output_dir = "streamlines_output"
os.makedirs(output_dir, exist_ok=True)

# Create subdirectories for each available pressure level
level_dirs = {level: os.path.join(output_dir, f"{int(level)} mb") for level in pressure_levels}
for dir_path in level_dirs.values():
    os.makedirs(dir_path, exist_ok=True)

# Define wind speed color levels and corresponding colors
cmap_colors = ["#dedede", "#95cee2", "#c1e4c4", "#bbdc71", "#8cce39", "#0eab42", "#fbee6b", "#eecf52", "#eea056", "#db4218"]
cmap_levels = [0, 1, 2, 4, 6, 9, 11, 14, 17, 21, 24]
cmap = mcolors.ListedColormap(cmap_colors)
norm = mcolors.BoundaryNorm(cmap_levels, cmap.N)

# Loop through time steps at 24-hour intervals
for time_step in dataset.valid_time[::24]:
    time_str = np.datetime_as_string(time_step.values, unit='h')  # Format time
    
    for pressure_level in pressure_levels:
        dataset_level = dataset.sel(pressure_level=pressure_level, valid_time=time_step,
                                    latitude=slice(lat_max, lat_min),  # Invert if needed
                                    longitude=slice(lon_min, lon_max))
        
        # Extract variables
        lon = dataset_level.longitude.values
        lat = dataset_level.latitude.values
        u_wind = dataset_level.u.values  # Selected time step
        v_wind = dataset_level.v.values  # Selected time step
        wind_speed = np.sqrt(u_wind**2 + v_wind**2)  # Compute wind speed
        
        # Create figure
        fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={"projection": ccrs.PlateCarree()})
        ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        
        # Add map features
        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.BORDERS, linestyle=':')
        ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray')
        gl = ax.gridlines(draw_labels=True, linestyle="--", linewidth=0.5)
        gl.top_labels = False  # Remove top latitude labels
        gl.right_labels = False  # Remove right longitude labels
        
        # Overlay wind speed
        wind_contour = ax.pcolormesh(lon, lat, wind_speed, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
        
        # Create streamplot with denser lines
        lon_grid, lat_grid = np.meshgrid(lon, lat)
        ax.streamplot(lon_grid, lat_grid, u_wind, v_wind, density=2.0, color="black", transform=ccrs.PlateCarree())
        
        # Add color bar with labels, making it 50% smaller and placing it on the right
        cbar = plt.colorbar(wind_contour, ax=ax, orientation="vertical", pad=0.02, aspect=25)
        cbar.set_label("Wind Speed (m/s)")
        cbar.set_ticks(cmap_levels)
        
        # Add title
        title = f"{int(pressure_level)} mb Streamlines - {time_str} UTC"
        ax.set_title(title)
        
        # Save plot
        output_filename = os.path.join(level_dirs[pressure_level], f"streamline_{int(pressure_level)}mb_{time_str}.png")
        plt.savefig(output_filename, dpi=300, bbox_inches='tight')
        plt.close(fig)
        
        # Print progress
        print(f"Finished outputting: {output_filename}")

# Close dataset
dataset.close()

# Final confirmation message
print("All streamline images successfully generated and saved.")