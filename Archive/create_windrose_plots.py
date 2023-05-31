from windrose import WindroseAxes
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import numpy as np
import os

import arcpy
arcpy.CheckOutExtension("GeoStats")
import datetime
from arcpy.sa import *

arcpy.env.overwriteOutput = True

arcpy.CheckOutExtension("Spatial")

wind_data_dir = r"\\loxodonta\gis\Source_Data\climatologyMeteorologyAtmosphere\state\CA\wind_data_ucsd"
#fire_polygons = r"E:\Projects\USGS-OAG_Wind_Temp_Drought_2018_mike_gough\Tasks\New_Wind_Data_Task\Data\Intermediate\TestData\TestData.gdb\tubbs_fire_project"
fire_polygons = r"E:\Projects\USGS-OAG_Wind_Temp_Drought_2018_mike_gough\Tasks\New_Wind_Data_Task\Data\Inputs\Inputs.gdb\firep19_1_project_GCS"
fire_polygons = r"E:\Projects\USGS-OAG_Wind_Temp_Drought_2018_mike_gough\Tasks\New_Wind_Data_Task_New_Data\Data\Intermediate\Examine_Wind_Data_Issues.gdb\Ranch_Fire_Largest_Fire_GCS_WGS84"

fire_date_field = "ALARM_DATE"
fire_name_field = "FIRE_NAME"
fire_area_field = "GIS_ACRES"
fire_area_units = "Acres"

plots_output_folder = r"E:\Projects\USGS-OAG_Wind_Temp_Drought_2018_mike_gough\Tasks\New_Wind_Data_Task\Tools\Scripts\Plots"

fire_year_start = 2018 
fire_year_end = 2019  # Will exclude fires that occurred on this year. 
oid_start = 0 

arcpy.env.workspace = r"in_memory"
#arcpy.env.workspace = r"c:\temp\temp.gdb"

# Won't include the last value. So 0, 1 will just be hour 0. NetCDF Hours just go up to 23, so set to 0-24 for the entire day. 
#12 am would be 0,1
#1 am would be 1,2
#2 am would be 2,3
#3 am would be 3,4
hour_chunks = [
    [0, 24],
]

directions = [
    [0, 360],
]


start_time = datetime.datetime.now()
print "\nStart Time: " + str(start_time)

fire_dates_processed = []

# Iterate over each fire polygon. Get the fire date and process fires for that date.
with arcpy.da.SearchCursor(fire_polygons, ["OBJECTID", fire_date_field, fire_name_field, fire_area_field]) as sc:

    for row in sc:
        object_id = row[0]
        fire_date = row[1]
        if not row[2]:
            fire_name = "UNNAMED"
        else:
            #fire with OID 20115 has the name SCHOOL\r\nSCHOOL. Screws up the output file path. 
            fire_name = row[2].split("\r\n")[0].replace(" ", "_")
        fire_area = row[3]

        wind_speed_array = []
        wind_dir_array = []

        # If the polygon has a valid fire date that is not None,
        # and we haven't already processed the rasters for this year
        # and the fire occurred after the user defined start year.
        if fire_date and fire_date not in fire_dates_processed and int(fire_date.year) >= fire_year_start and int(fire_date.year) < fire_year_end and object_id >= oid_start:

            query ='"OBJECTID" = ' + str(object_id)
            arcpy.MakeFeatureLayer_management(fire_polygons, "fire_polygon", query) 

            # Add the fire date to the list of fire dates already processed (rasters have already been created).
            fire_dates_processed.append(fire_date)

            fire_year = (str(fire_date.year)).zfill(2)
            fire_month = (str(fire_date.month)).zfill(2)
            fire_day = (str(fire_date.day)).zfill(2)

            print "Fire OBJECTID: " + str(object_id)
            print "Fire Name: " + fire_name
            print "Fire Date: " + str(fire_date).split(" ")[0]
            print "Fire Area: " + str(fire_area) + " " + fire_area_units

            # Get the path to the netCDF for the year of the fire.
            netcdf_file_name = "LOCA-ERA5.2020-04-22.wind_spd_dir.%s.pdbc.nc" % fire_year
            input_netcdf_file = wind_data_dir + os.sep + netcdf_file_name

            # Iterate over each of the hour chunks. Typically just one. 
            for hour_chunk in hour_chunks:

                start_hour = hour_chunk[0]
                end_hour = hour_chunk[1]

                print "\nCalculating Hour Chunk"
                print "Start Hour: " + str(start_hour)
                print "End Hour: " + str(end_hour)

                raster_dict = {}

                # Create a wind direction raster layer and a wind speed raster layer for each hour between the start hour and the end hour for this chunk.
                #print "\nCreating NetCDF Raster Layers (One for Wind Direction, One for Wind Speed)"
                for hour in range(start_hour, end_hour):
                    
                    print "Hour: " + str(hour)

                    time = '%s/%s/%s %s:00:00' % (fire_month, fire_day, fire_year, hour)
                    print "Time Dimension: " + str(time)
                    dimension_value = "time '%s'" % time

                    # Wind Direction NetCDF Layer
                    netcdf_raster_layer_wind_dir = "wind_dir_" + "_".join([str(fire_year), str(fire_month), str(fire_day), str(hour)])
                    arcpy.MakeNetCDFRasterLayer_md(
                        in_netCDF_file=input_netcdf_file,
                        variable="wind_dir", x_dimension="lon", y_dimension="lat", out_raster_layer=netcdf_raster_layer_wind_dir,
                        band_dimension="", dimension_values=dimension_value, value_selection_method="BY_VALUE")

                    # Wind Speed NetCDF Layer
                    netcdf_raster_layer_wind_speed = "wind_speed_" + "_".join([str(fire_year), str(fire_month), str(fire_day), str(hour)])
                    arcpy.MakeNetCDFRasterLayer_md(
                        in_netCDF_file=input_netcdf_file,
                        variable="wind_speed", x_dimension="lon", y_dimension="lat", out_raster_layer=netcdf_raster_layer_wind_speed,
                        band_dimension="", dimension_values=dimension_value, value_selection_method="BY_VALUE")

                    # Get wind speed pixels for north facing pixels for this hou. Add raster path to the list.
                    #wind_speed_list = arcpy.sa.Con(netcdf_raster_layer_wind_dir, netcdf_raster_layer_wind_speed, 0, "VALUE >= 0")
                    #wind_speed_list.append(wind_speed_n)
                    print "Extracting Wind Speed Values..."
                    arcpy.ExtractValuesToTable_ga("fire_polygon", netcdf_raster_layer_wind_speed, "temp_table_wind_speed")
                    with arcpy.da.SearchCursor("temp_table_wind_speed", ["Value"]) as sc2:
                        for row2 in sc2:
                            print row2[0]
                            wind_speed_array.append(row2[0])
                            
                        # Check to make sure there's wind speed data. If not not, break out of the for loop. 
                        if len(wind_speed_array) == 0:
                            print "No wind speed data in the search cursor. Skipping this fire..."
                            break
                            
                    print "Extracting Wind Direction Values..."
                    arcpy.ExtractValuesToTable_ga("fire_polygon", netcdf_raster_layer_wind_dir, "temp_table_wind_dir")
                    with arcpy.da.SearchCursor("temp_table_wind_dir", ["Value"]) as sc3:
                        for row3 in sc3:
                            print row3[0]
                            wind_dir_array.append(row3[0])


            if len(wind_speed_array)> 0 and len(wind_dir_array) > 0:
                print "\nCreating windrose plot..."
                # Create windrose plots. 
                #ws = np.random.random(500) * 6
                #wd = np.random.random(500) * 360
                ws = np.array(wind_speed_array)
                wd = np.array(wind_dir_array)
                ax = WindroseAxes.from_ax()
                # It looks like the default angular labels are wrong!
                #ax.set_xticklabels(['N', 'NW', 'W', 'SW', 'S', 'SE', 'E', 'NE'])
                # Labels are off in windrose: https://github.com/python-windrose/windrose/issues/151 
                #ax.set_xticklabels([90, 45, 0, 315, 270, 225, 180, 135])
                ax.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'])
                ax.bar(wd, ws, normed=True, opening=0.8, edgecolor='white')
                ax.set_legend(title='Wind Speed in m/s')
                    
                plot_title = fire_name + " (" + str(fire_date).split(" ")[0] + ", " + str(int(fire_area)) + " " + fire_area_units + ")"     
                plt.title(plot_title)
                
                output_file = plots_output_folder + os.sep + fire_name + "_" + str(fire_date).split(" ")[0] + "_" + str(object_id) + ".png"
                plt.savefig(output_file)

end_time = datetime.datetime.now()
print "\nEnd Time: " + str(end_time)
print "Duration: " + str(end_time - start_time)