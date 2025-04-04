import os
import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import geopandas as gpd
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def find_nc_files(folder):
    """Retrieve all .nc files in a given folder."""
    return [os.path.join(folder, f) for f in os.listdir(folder) if f.endswith(".nc")]

# Ensure output directory exists
output_dir = "bias_gsmap-station"
os.makedirs(output_dir, exist_ok=True)

# Load station rainfall data
station_file = "rainfall_data.csv"
df_station = pd.read_csv(station_file)
df_station.columns = df_station.columns.str.strip()

# Identify rainfall columns dynamically (columns with YYYYMMDD format)
date_columns = [col for col in df_station.columns if col.isdigit() and len(col) == 8]

df_station["lat"] = pd.to_numeric(df_station["lat"], errors="coerce")
df_station["lon"] = pd.to_numeric(df_station["lon"], errors="coerce")

# Load the shapefile
shapefile = gpd.read_file("shapefiles/phprov.shp")

# Define bias color scale
levels = [-100, -50, -25, -10, -5, -1, 0, 1, 5, 10, 25, 50, 100]
colors = ["#a80000", "#ff0000", "#ff7373", "#ffb8b8", "#ffe6e6", "#ffffff", "#e6f7ff", "#b8e6ff", "#73c2ff", "#7373ff", "#0000a8", "#000066"]
cmap = mcolors.ListedColormap(colors)
norm = mcolors.BoundaryNorm(levels, cmap.N)

# Iterate through each date column
for date in date_columns:
    station_data = df_station[["lat", "lon", date]].copy()
    station_data[date] = station_data[date].replace("T", 0).astype(float)
    
    # Path to the GSMaP folder for this date
    gsmap_folder = os.path.join("gsmap_data", date)
    if not os.path.exists(gsmap_folder):
        print(f"Skipping {date}, GSMaP data folder not found.")
        continue
    
    # Accumulate 24-hour rainfall from all available NetCDF files
    gsmap_total = None
    nc_files = find_nc_files(gsmap_folder)
    if not nc_files:
        print(f"Skipping {date}, no GSMaP NetCDF files found.")
        continue
    
    for nc_file in nc_files:
        ds = xr.open_dataset(nc_file, decode_times=False)
        
        # Extract correct rainfall variable
        if "hourlyPrecipRate" in ds:
            rainfall = ds["hourlyPrecipRate"].values.squeeze()
        elif "hourlyPrecipRateGC" in ds:
            rainfall = ds["hourlyPrecipRateGC"].values.squeeze()
        else:
            print(f"Skipping {nc_file}, no valid precipitation variable found.")
            ds.close()
            continue
        
        lats = ds["Latitude"].values
        lons = ds["Longitude"].values
        
        if gsmap_total is None:
            gsmap_total = rainfall
        else:
            gsmap_total += rainfall
        
        ds.close()
    
    if gsmap_total is None:
        print(f"Skipping {date}, no valid GSMaP data found.")
        continue
    
    # Compute bias: GSMaP 24-hour rainfall - Station rainfall
    station_data["bias"] = np.nan
    for i, row in station_data.iterrows():
        lat_idx = np.abs(lats - row["lat"]).argmin()
        lon_idx = np.abs(lons - row["lon"]).argmin()
        station_data.at[i, "bias"] = gsmap_total[lat_idx, lon_idx] - row[date]
    
    # Plot
    fig, ax = plt.subplots(figsize=(8, 6), subplot_kw={"projection": ccrs.PlateCarree()})
    ax.set_extent([df_station["lon"].min() - 1, df_station["lon"].max() + 1, df_station["lat"].min() - 1, df_station["lat"].max() + 1], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.BORDERS, linestyle=":")
    ax.add_feature(cfeature.COASTLINE)
    shapefile.boundary.plot(ax=ax, edgecolor="black", linewidth=1, transform=ccrs.PlateCarree())
    
    # Scatter plot
    sc = ax.scatter(
        station_data["lon"], station_data["lat"], 
        c=station_data["bias"], cmap=cmap, norm=norm, 
        edgecolor="black", s=50, transform=ccrs.PlateCarree()
    )
    
    # Add color bar
    cbar = plt.colorbar(sc, ax=ax, orientation="vertical", shrink=0.7, pad=0.02, extend="both")
    cbar.set_label("Bias (GSMaP - Station) mm")
    cbar.set_ticks(levels)
    cbar.set_ticklabels([str(l) for l in levels])
    
    # Format date
    formatted_date = pd.to_datetime(date, format="%Y%m%d").strftime("%Y-%m-%d")
    ax.set_title(f"GSMaP - Station Rainfall Bias ({formatted_date})")
    
    # Longitude and latitude bars
    ax.set_xticks(np.arange(df_station["lon"].min() - 1, df_station["lon"].max() + 2, 5), crs=ccrs.PlateCarree())
    ax.set_yticks(np.arange(df_station["lat"].min() - 1, df_station["lat"].max() + 2, 5), crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, pos: f"{v:.1f}°E"))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, pos: f"{v:.1f}°N"))
    ax.tick_params(axis="both", labelsize=10)
    ax.set_xlabel("")
    ax.set_ylabel("")
    
    # Save
    plt.savefig(f"{output_dir}/{date}.png", dpi=300, bbox_inches="tight")
    plt.close()
    
print(f"Bias plots saved in '{output_dir}'")
