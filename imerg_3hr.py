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

# Group into 3-hour chunks (6 files)
group_size = 6
results = []

for i in range(0, len(nc_files), group_size):
    group = nc_files[i:i + group_size]
    if len(group) < group_size:
        break  # Skip incomplete groups

    rainfall_total = 0.0
    start_time = None

    for file in group:
        try:
            ds = xr.open_dataset(file)
            rain_data = ds['precipitation']
            rainfall = rain_data.sel(lat=target_lat, lon=target_lon, method='nearest').values.item()
            rainfall_total += rainfall
            if start_time is None:
                start_time = np.datetime64(ds['time'].values[0])
            ds.close()
        except Exception as e:
            print(f"Error processing {file}: {e}")
            break

    if start_time is not None:
        ts_str = str(start_time)
        mmddhhmm = f"{ts_str[5:7]}-{ts_str[8:10]}-{ts_str[11:13]}{ts_str[14:16]}"
        results.append([mmddhhmm, f"{rainfall_total:.2f}"])

# Write to CSV
csv_file = 'imerg_3hr.csv'
with open(csv_file, mode='w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Timestamp', 'Rainfall_mm'])
    writer.writerows(results)

print(f"Saved rainfall data to {csv_file}")
