# Print dataset dimensions before slicing
print(f"Processing {nc_file}")
print("Before slicing:", dataset.dims)

# Select the specified domain first
dataset = dataset.sel(Latitude=slice(lat_max, lat_min), Longitude=slice(lon_min, lon_max))

# Print dataset dimensions after slicing
print("After slicing:", dataset.dims)
