from netCDF4 import Dataset

# Open the NetCDF file
nc_file = "era5/8c426675800d082614e8c529e1d8b0c1.nc"
data = Dataset(nc_file)

# Print the names and dimensions of all variables to check for time variable
for var_name in data.variables:
    print(f"Variable: {var_name}")
    print(f"Dimensions: {data.variables[var_name].dimensions}")
    print(f"Shape: {data.variables[var_name].shape}")
    print(f"Attributes: {data.variables[var_name].__dict__}")
    print()

# Close the dataset
data.close()
