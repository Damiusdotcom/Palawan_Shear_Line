import xarray as xr

# Set the path to your NetCDF file
file_path = 'era5/8c426675800d082614e8c529e1d8b0c1.nc'  # Replace with the actual path

# Open the dataset
ds = xr.open_dataset(file_path)

# Print basic dataset info
print("\n📦 Dataset Info:\n")
print(ds)

# Print variable names
print("\n📌 Variables in the dataset:\n")
for var in ds.variables:
    print(f" - {var}")

# Optional: Show attributes for each variable
print("\n📋 Variable Attributes:\n")
for var in ds.variables:
    print(f"\n🔹 {var}")
    print(ds[var].attrs)

# Close the dataset
ds.close()
