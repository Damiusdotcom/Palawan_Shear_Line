import xarray as xr

# Set the path to your NetCDF file
#file_path = 'Y:\\wfs_shared\\Personnel Files\\DJEV\\Palawan_Shear_Line\\mean_streamlines\\v_2020.nc'  # Replace with the actual path

file_path = 'era5\\era5_file.nc'  # Replace with the actual path

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
