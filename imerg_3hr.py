import os
import xarray as xr
import numpy as np
import csv

# Folder containing 30-minute IMERG NetCDF files
data_folder = 'imerg_30min'

# Coordinate of interest
target_lat = 9.740134
target_lon = 118.758613

# List and sort NetCDF files by internal time variable
def get_file_time(filepath):
    try:
        with xr.open_dataset(filepath) as ds:
            return np.datetime64(ds['time'].values[0])
    except Exception as e:
        print(f"Error reading time from {filepath}: {e}")
        return None

nc_files = sorted(
    [os.path.join(data_folder, f) for f in os.listdir(data_folder) if f.endswith('.nc4')],
    key=get_file_time
)

# Prepare to collect results
results = []

# Process each file individually
for file in nc_files:
    try:
        ds = xr.open_dataset(file)
        rain_data = ds['precipitation']
        rainfall = rain_data.sel(lat=target_lat, lon=target_lon, method='nearest').values.item()
        timestamp = str(ds['time'].values[0])
        print(f"{timestamp}: {rainfall:.2f} mm")
        results.append([timestamp, f"{rainfall:.2f}"])
        ds.close()
    except Exception as e:
        print(f"Error processing {file}: {e}")

# Write to CSV
csv_file = 'imerg_30min.csv'
with open(csv_file, mode='w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Timestamp', 'Rainfall_mm'])
    writer.writerows(results)

print(f"Saved rainfall data to {csv_file}")
