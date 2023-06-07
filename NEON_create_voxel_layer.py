import os
import arcpy
import netCDF4
import numpy as np
import datetime as dt
#import arcpy
#input_gridded_points = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Test\Test.gdb\lidar_multi_points_to_points_join_xyz_xy_table_to_point"
#output_gridded_points = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Test\Test.gdb\lidar_multi_points_to_points_join_xyz_xy_table_to_point_output"

#arcpy.CopyFeatures_management(input_gridded_points, output_gridded_points)

#with arcpy.UpdateCursor(output_gridded_points, ["Z_MAX"]) as uc:
#    for row in uc:
#        print(row[0])

filePath = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Volume\Points_To_CSV\lidar_multi_points_to_points_join_xyz_xy_table_to_point_subset.csv"
output_netcdf = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Volume\NetCDF\csv_to_netcdf.nc"

from netCDF4 import Dataset
import numpy as np
import pandas as pd

#Create a pandas dataframe and insert data from CSV/TEXT file
dfPoints = pd.read_csv(filePath)

#Sort values to ensure they are in the correct order
#dfPoints = dfPoints.sort_values(by=['T','Z','Y','X'])
dfPoints = dfPoints.sort_values(by=['Date', 'Z_Max','POINT_Y','POINT_X'])

#Create domain for longitude/latitude
#Each domain has unique values, no repeating numbers, and are sorted (to be monotonic)
xDomain = np.sort(np.unique(dfPoints.iloc[:,8].values)) # 0th column contains x values
print(xDomain)
yDomain = np.sort(np.unique(dfPoints.iloc[:,9].values)) # 1st column contains y values
print(yDomain)
zDomain = np.sort(np.unique(dfPoints.iloc[:,12].values)) # 2nd column contains z values
print(zDomain)
tDomain = np.sort(np.unique(dfPoints.iloc[:,10].values)) # 3rd column contains t values
print(tDomain)

#Create NetCDF file
outDataSet = Dataset(output_netcdf, 'w', format = 'NETCDF4') # Creates the output NetCDF file
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
    #dfPoints['Data'].values,
    test,
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
#outDataSet.esri_pe_string = 'PROJCS["ETRS_1989_UTM_Zone_32N_7stellen",GEOGCS["GCS_ETRS_1989",DATUM["D_ETRS_1989",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",2500000.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",9.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0],AUTHORITY["Esri",102328]]'
outDataSet.esri_pe_string = 'PROJCS["WGS_1984_UTM_Zone_11N",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",-117.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]],VERTCS["unknown",VDATUM["unknown"],PARAMETER["Vertical_Shift",0.0],PARAMETER["Direction",1.0],UNIT["Meter",1.0]]',

outDataSet.close()

