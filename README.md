The case study will analyze the rainfall patterns of the heavy rainfall event in Puerto Princesa Palawan on February 7-11, 2025.

Plots are overlayed over shapefiles/phprov.shp

Domain is 4N to 22N, 114E to 130E

Directories and Files:

imerg/ <br>
    raw netcdf file/s of IMERG

era5/<br>
    raw nc file/s of era5

imerg_plots/ <br>
    output plots of IMERG variable "precipitation"
    
imerg_point/<br>
    point plots of IMERG variable "precipitation" 
    
rainfall_actual/<br>
    point plots of actual recorded rainfall by PAGASA synoptic stations
    
bias_imerg-station/<br>
    bias plot, formula is imerg - station rainfall

streamlines_output/
    streamline output using the nc file inside era5/ as input
    
coordinates.csv<br>
    latitude and longitude of PAGASA synoptic stations. Column titles should be 'lat' and 'lon'
    
rainfall_data.csv<br>
    Extracted rainfall data from PAGASA-WD Sycoder. Column titles should be yyyymmdd, 'lat', and 'lon'
    
imerg_plot.py, imerg_point.py, decoded.py, imerg_station.py<br>
    Python scripts to output rainfall charts<br>
    The scripts will automatically use all netcdf files under the imerg directory as input<br>
    The scripts will also use rainfall data under the rainfall_data.csv as input<br>
<br>

streamlines.py<br>
    Python script to output streamlines using the nc file inside era5/<br>
    Subfolders for each available level will be automatically created under streamlines/output if no subfolders are detected <br>
    Note: Calling the file is still not automatic, nc_file variable should be changed manually <br>

For this case, 5 days of rainfall data was used, but the script can manage more data as long as the netcdf files are placed in the imerg directory and the column title format in rainfal_data.csv is correct.
 