The case study will analyze the rainfall patterns of the heavy rainfall event in Puerto Princesa Palawan on February 7-11, 2025.

Plots are overlayed over shapefiles/phprov.shp

Domain is 4N to 22N, 114E to 130E

Directories and Files:

imerg <br>
    raw netcdf files of IMERG
    
imerg_plots<br>
    output plots of IMERG variable "precipitation"
    
imerg_point<br>
    point plots of IMERG variable "precipitation" 
    
rainfall_actual<br>
    point plots of actual recorded rainfall by PAGASA synoptic stations
    
bias_imerg-station<br>
    bias plot, formula is imerg - station rainfall
    
coordinates.csv<br>
    latitude and longitude of PAGASA synoptic stations. Column titles should be 'lat' and 'lon'
    
rainfall_data.csv<br>
    Extracted rainfall data from PAGASA-WD Sycoder. Column titles should be yyyymmdd, 'lat', and 'lon'
    
imerg_plot.py, imerg_point.py, decoded.py, imerg_station.py<br>
    Python scripts to output the charts<br>
    The scripts will automatically use all netcdf files under the imerg directory as input<br>
    The scripts will also use rainfall data under the rainfall_data.csv as input<br>
<br>
For this case, 5 days of rainfall data was used, but the script can manage more data as long as the netcdf files are placed in the imerg directory and the column title format in rainfal_data.csv is correct.
 