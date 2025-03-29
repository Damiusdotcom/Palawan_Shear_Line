# script created by DJEV
import os
import ssl
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import geopandas as gpd
import cartopy.crs as ccrs
import cartopy.feature as cfeature

ssl._create_default_https_context = ssl._create_unverified_context

# Ensure output directory exists
output_dir = "rainfall_actual"
os.makedirs(output_dir, exist_ok=True)

# Load CSV file
csv_file = "rainfall_data.csv"
df = pd.read_csv(csv_file)

# Ensure column names are stripped of whitespace
df.columns = df.columns.str.strip()

# Identify rainfall columns dynamically (columns with YYYYMMDD format)
date_columns = [col for col in df.columns if col.isdigit() and len(col) == 8]

# Convert lat/lon to numeric
df["lat"] = pd.to_numeric(df["lat"], errors="coerce")
df["lon"] = pd.to_numeric(df["lon"], errors="coerce")

# Load the shapefile
shapefile = gpd.read_file("shapefiles/phprov.shp")

# Define precipitation thresholds (mm) and corresponding colors
levels = [0, 1, 10, 25, 50, 100, 200, 300, 500, 700]
colors = ["#ffffff", "#bab8b8", "#00c5ff", "#6bfb90", "#ffff00", "#ffaa00", "#ff0000", "#ff73df", "#8400a8"]

# Create a custom colormap and norm
cmap = mcolors.ListedColormap(colors).with_extremes(over="#000000")
norm = mcolors.BoundaryNorm(levels, cmap.N)

# Loop through each date column to create a plot
for date in date_columns:
    # Replace "T" (trace amounts) with 0 and convert to numeric
    df[date] = df[date].replace("T", 0)
    df[date] = pd.to_numeric(df[date], errors="coerce")

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6), subplot_kw={"projection": ccrs.PlateCarree()})
    
    # Set map extent based on station locations
    ax.set_extent([df["lon"].min() - 1, df["lon"].max() + 1, df["lat"].min() - 1, df["lat"].max() + 1], crs=ccrs.PlateCarree())

    # Add map features
    ax.add_feature(cfeature.BORDERS, linestyle=":")
    ax.add_feature(cfeature.COASTLINE)

    # Plot shapefile
    shapefile.boundary.plot(ax=ax, edgecolor="black", linewidth=1, transform=ccrs.PlateCarree())

    # Scatter plot with precipitation data
    sc = ax.scatter(
        df["lon"], df["lat"],  
        c=df[date], cmap=cmap, norm=norm,
        edgecolor="black", s=50, transform=ccrs.PlateCarree()
    )

    # Add color bar
    cbar = plt.colorbar(sc, ax=ax, orientation="vertical", shrink=0.7, pad=0.02, extend="max")
    cbar.set_label("Precipitation (mm)")
    cbar.set_ticks(levels)
    cbar.set_ticklabels([str(l) for l in levels])

    # Format date
    formatted_date = pd.to_datetime(date, format="%Y%m%d").strftime("%Y-%m-%d")

    # Set title
    ax.set_title(f"Station 24-hour Accumulated Rainfall ({formatted_date})", fontsize=14)

    # Define tick values for longitude and latitude at 5-degree increments
    lon_ticks = np.arange(df["lon"].min() - 1, df["lon"].max() + 2, 5)
    lat_ticks = np.arange(df["lat"].min() - 1, df["lat"].max() + 2, 5)

    # Set tick positions
    ax.set_xticks(lon_ticks, crs=ccrs.PlateCarree())
    ax.set_yticks(lat_ticks, crs=ccrs.PlateCarree())

    # Format tick labels to display degrees
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, pos: f"{v:.1f}°E"))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, pos: f"{v:.1f}°N"))

    # Adjust tick label size
    ax.tick_params(axis="both", labelsize=10)
    
    # Remove labels
    ax.set_xlabel("")
    ax.set_ylabel("")

    # Save the plot
    plt.savefig(f"{output_dir}/{date}.png", dpi=300, bbox_inches="tight")
    plt.close()

print(f"Plots saved in '{output_dir}'")
