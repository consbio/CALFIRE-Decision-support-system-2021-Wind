########################################################################################################################
# File name: NEON_create_voxel_layer.py
# Author: Mike Gough
# Date created: 06/06/2023
# Python Version: 2.7
# Description:
# Creates a NetCDF Voxel Layer from a CSV file generated from a vector fishnet containing point counts at different
# heights above the ground.
########################################################################################################################
import os
import arcpy
from datetime import datetime
import math
import netCDF4
import numpy as np
import datetime as dt
arcpy.env.overwriteOutput = True

tmp_gdb = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Scratch\Scratch.gdb"

output_csv = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Volume\CSV\merge_fc.csv"

output_netcdf_voxel = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\NetCDF_Voxel\csv_to_netcdf_voxel.nc"

#count_lidar_point_returns_output = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\Outputs.gdb\segmented_heights_above_ground_Extent_1_Tower_Location_subset"
count_lidar_point_returns_output = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\Outputs.gdb\segmented_heights_above_ground_Extent_1_Tower_Location"

start_script = datetime.now()
print("\nStart Time: " + str(start_script))

def create_csv():

    elevation_slice_fields = [field.name for field in arcpy.ListFields(count_lidar_point_returns_output) if "count_" in field.name]

    # Convert each field slice to a separate feature class
    print("Converting each field slice to a separate feature class...")
    merge_list = []
    for e_field in elevation_slice_fields:
        print(e_field)
        tmp_conversion_fc = os.path.join(tmp_gdb, e_field + "_tmp")
        #fields_to_export = "'" + field +  '"' + field + '"' + "true true false 4 long"

        # Create a field map of desired output fields.
        myfields = [e_field, 'GRID_ID', 'DTM_Extraction']
        # create an empty field mapping object
        mapS = arcpy.FieldMappings()
        #for each field, create an individual field map, and add it to the field mapping object
        for field in myfields:
            map = arcpy.FieldMap()
            map.addInputField(count_lidar_point_returns_output, field)
            mapS.addFieldMap(map)

        arcpy.conversion.ExportFeatures(count_lidar_point_returns_output, tmp_conversion_fc, '', "NOT_USE_ALIAS", mapS, None)

        elevation_slice_min = e_field.split("_")[-1].replace("m", "")
        arcpy.AddField_management(tmp_conversion_fc, "Z_Slice", "LONG")
        arcpy.CalculateField_management(tmp_conversion_fc, "Z_Slice", elevation_slice_min)
        arcpy.AlterField_management(tmp_conversion_fc, e_field, "count", "count")

        merge_list.append(tmp_conversion_fc)

    print("Merging...")
    merge_fc = os.path.join(tmp_gdb, "tmp_merge")
    arcpy.Merge_management(merge_list, merge_fc)

    print("Adding an X field...")
    arcpy.AddField_management(merge_fc, "X", "Double")

    print("Adding a Y field...")
    arcpy.AddField_management(merge_fc, "Y", "Double")

    print("Calculating X and Y Centroid Coordinates...")
    arcpy.management.CalculateGeometryAttributes(
        merge_fc,
        "X CENTROID_X;Y CENTROID_Y", '', '', None, "SAME_AS_INPUT")

    print("Adding a time field...")
    # Add a time field. Pseudo-time of 1 for now.
    arcpy.AddField_management(merge_fc, "Time", "SHORT")
    arcpy.CalculateField_management(merge_fc, "Time", 1)

    print("Creating the output CSV..")
    output_csv_dir = os.path.dirname(output_csv)
    output_csv_name = os.path.basename(output_csv)
    arcpy.TableToTable_conversion(merge_fc, output_csv_dir, output_csv_name)

def create_voxel():

    #filePath = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Volume\Points_To_CSV\lidar_multi_points_to_points_join_xyz_xy_table_to_point_subset.csv"

    from netCDF4 import Dataset
    import numpy as np
    import pandas as pd

    #Create a pandas dataframe and insert data from CSV/TEXT file
    dfPoints = pd.read_csv(output_csv)

    #Sort values to ensure they are in the correct order
    #dfPoints = dfPoints.sort_values(by=['T','Z','Y','X'])
    #dfPoints = dfPoints.sort_values(by=['Date', 'Z_Max','POINT_Y','POINT_X'])
    dfPoints = dfPoints.sort_values(by=['Time', 'Z_Slice','Y','X'])

    #Create domain for longitude/latitude
    #Each domain has unique values, no repeating numbers, and are sorted (to be monotonic)
    xDomain = np.sort(np.unique(dfPoints.iloc[:,8].values)) # 0th column contains x values
    print(xDomain)
    yDomain = np.sort(np.unique(dfPoints.iloc[:,9].values)) # 1st column contains y values
    print(yDomain)
    zDomain = np.sort(np.unique(dfPoints.iloc[:,4].values)) # 2nd column contains z values
    print(zDomain)
    tDomain = np.sort(np.unique(dfPoints.iloc[:,7].values)) # 3rd column contains t values
    print(tDomain)

    #Create NetCDF file
    outDataSet = Dataset(output_netcdf_voxel, 'w', format = 'NETCDF4') # Creates the output NetCDF file
    outDataSet.createDimension('x', len(xDomain))                # Creates the x dimension
    outDataSet.createDimension('y', len(yDomain))                # Creates the y dimension
    outDataSet.createDimension('z', len(zDomain))                # Creates the z dimension
    outDataSet.createDimension('t', len(tDomain))                # Creates the t dimension

    #Create variables
    ncX = outDataSet.createVariable('x', np.float32, 'x') # Create variable x
    ncY = outDataSet.createVariable('y', np.float32, 'y') # Create variable y
    ncZ = outDataSet.createVariable('z', np.float32, 'z') # Create variable z
    ncT = outDataSet.createVariable('t', np.float32, 't') # Create variable t
    #ncT = outDataSet.createVariable('t', np.datetime64, 't') # Create variable t

    #Create variable data with dimensions (t,z,y,x). The fill value is set to -99999 which are values to be ignored by client.
    ncData = outDataSet.createVariable('data',np.float32,('t','z','y','x'),fill_value = -99999)

    #Assign values1
    ncX[:] = xDomain[:]
    ncY[:] = yDomain[:]
    ncZ[:] = zDomain[:]
    ncT[:] = tDomain[:]

    test = np.arange(900)
    #should be 900 points
    #The dataframe 'Data' column must be reshaped to match the dimension shapes and placed into the ncData variable
    ncData[:,:,:,:]  = np.reshape(
        dfPoints['count'].values,
        #test,
        (tDomain.shape[0],zDomain.shape[0],yDomain.shape[0],xDomain.shape[0])
        )

    #Assign variable attributes

    ncData.long_name = "My test data"

    ncZ.positive = 'up'

    ncX.standard_name = 'projection_x_coordinate'
    ncX.units = 'm'

    ncY.standard_name = 'projection_y_coordinate'
    ncY.units = 'm'

    #ncT.units = 'hours since 2000-01-01 00:00:00'

    #Assign global attribute. This attribute is to assign a coordinate system.
    outDataSet.esri_pe_string = 'PROJCS["WGS_1984_UTM_Zone_11N",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",-117.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]],VERTCS["unknown",VDATUM["unknown"],PARAMETER["Vertical_Shift",0.0],PARAMETER["Direction",1.0],UNIT["Meter",1.0]]',

    outDataSet.close()

create_csv()
create_voxel()

end_script = datetime.now()
print("\nEnd Time: " + str(end_script))
print("Duration: " + str(end_script - start_script))


