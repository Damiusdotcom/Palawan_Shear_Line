import xarray as xr

# File path to the NetCDF file
file_path = 'era5/8c426675800d082614e8c529e1d8b0c1.nc'

# Load the NetCDF file using xarray
ds = xr.open_dataset(file_path)

# Print all variables and their dimensions
print("Variables in the NetCDF file:")
for var in ds.variables:
    print(f"Variable: {var}")
    print(f"Dimensions: {ds[var].dims}")
    print(f"Shape: {ds[var].shape}")
    print(f"Attributes: {ds[var].attrs}\n")

# Optionally, print the dimensions of the entire dataset
print("\nDimensions of the dataset:")
for dim in ds.dims:
    print(f"Dimension: {dim}, Size: {ds.dims[dim]}")

# Close the dataset
ds.close()
