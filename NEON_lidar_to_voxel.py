########################################################################################################################
# File name: NEON_lidar_to_voxel.py
# Author: Mike Gough
# Date created: 06/13/2023
# Python Version: 3.x
# Description:
# Creates a NetCDF voxel layer from a LiDAR file and an extent. The voxel layer contains an aggregation of LiDAR point
# counts on a regular volumetric grid.
########################################################################################################################

import os
import arcpy
from datetime import datetime
from netCDF4 import Dataset
import numpy as np
import pandas as pd

arcpy.env.overwriteOutput = True

# Input Parameters

version = "v1"
#extent_fc = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Inputs\Extents\Extents.gdb\extent_1_tower_location"
extent_fc = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Inputs\Extents\Extents.gdb\SOAP_tile_298000_4100000"

extent_name = extent_fc.split(os.sep)[-1]

voxel_size = 1  # Used to define the x,y, and y dimensions of the voxel (units are m).
max_chm_offset = 10  # Max Canopy Height Offset allowed. LiDAR point returns greater than this distance above the CHM will be considered errors and will be deleted (units are m).
starting_height = 1  # Used to remove ground points. Any points less than this distance from the ground will not be included (units are m).

output_proj = 'PROJCS["WGS_1984_UTM_Zone_11N",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",-117.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]],VERTCS["unknown",VDATUM["unknown"],PARAMETER["Vertical_Shift",0.0],PARAMETER["Direction",1.0],UNIT["Meter",1.0]]'

NEON_lidar_laz_file = r"\\loxodonta\gis\Source_Data\environment\region\NEON_SITES\SOAP\Discrete_return_LiDAR_point_cloud\2019\06\NEON_lidar-point-cloud-line\NEON.D17.SOAP.DP1.30003.001.2019-06.basic.20230523T232633Z.RELEASE-2023\NEON_D17_SOAP_DP1_298000_4100000_classified_point_cloud_colorized.laz"
NEON_DTM = r"\\loxodonta\gis\Source_Data\environment\region\NEON_SITES\SOAP\Elevation_LiDAR\2021\07\NEON_lidar-elev\NEON.D17.SOAP.DP3.30024.001.2021-07.basic.20230601T181117Z.RELEASE-2023\NEON_D17_SOAP_DP3_298000_4100000_DTM.tif"
NEON_CHM = r"\\loxodonta\gis\Source_Data\environment\region\NEON_SITES\SOAP\Ecosystem_structure\2019\06\NEON_struct-ecosystem\NEON.D17.SOAP.DP3.30015.001.2019-06.basic.20230524T172838Z.RELEASE-2023\NEON_D17_SOAP_DP3_298000_4100000_CHM.tif"

tmp_gdb = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Scratch\Scratch.gdb"
tmp_folder = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Scratch"
csv_folder = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Volume\CSV"
intermediate_gdb = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Volume\Volume.gdb"
intermediate_folder = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Volume"

# Derived Parameters

lidar_folder = os.path.join(intermediate_folder, "Lidar")
lidar_file_conversion_name = NEON_lidar_laz_file.split(os.sep)[-1].split(".")[0] + "_lidar_las_file_conversion.lasd"
lidar_file_conversion = os.path.join(lidar_folder, lidar_file_conversion_name)

lidar_clip_file = os.path.join(lidar_folder, lidar_file_conversion_name.split(".")[0] + "_clip" + ".lasd")
input_las_file = os.path.join(lidar_folder, NEON_lidar_laz_file.split(os.sep)[-1].replace(".laz", ".las"))

lidar_multipoint = os.path.join(intermediate_gdb, "lidar_multipoint_" + extent_name + "_" + version)
lidar_multipoint_to_points = os.path.join(intermediate_gdb, "lidar_multipoint_to_points_" + extent_name + "_" + version)
output_dtm = os.path.join(intermediate_gdb, "neon_dtm_clipped_" + extent_name + "_" + version)
output_chm = os.path.join(intermediate_gdb, "neon_chm_clipped_" + extent_name + "_" + version)
lidar_points_with_dtm_extraction = os.path.join(tmp_gdb, "lidar_points_with_dtm_extraction_" + extent_name + "_" + version)
lidar_points_with_z_dtm_and_chm = os.path.join(intermediate_gdb, "lidar_points_with_z_dtm_and_chm_" + extent_name + "_" + version)

fishnet_input_points_extent = os.path.join(intermediate_gdb, "fishnet_lidar_point_" + extent_name + "_" + version)
merged_fishnet_points = os.path.join(intermediate_gdb, "merged_fishnet_points_" + extent_name + "_" + version)

# Output Files

output_csv = csv_folder + os.sep + "merge_fc_" + extent_name + "_" + version + ".csv"
output_netcdf_voxel = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\NetCDF_Voxel\csv_to_netcdf_voxel_ " + extent_name + "_" + version + ".nc"


########################################################################################################################

start_script = datetime.now()
print("\nStart Time: " + str(start_script))

las_file_parse_date = datetime.strptime(NEON_lidar_laz_file.split("NEON.")[-1].split(".")[5] + '+0000', "%Y-%m%z") # Specify Time Zone as GMT

print("\nLiDAR File Datetime: " + str(las_file_parse_date))

# String time
#las_file_date = las_file_parse_date.strftime("%Y-%m-%d %H:%M:%S")

# Epoch time
las_file_date = las_file_parse_date.timestamp()

print("LiDAR File Epoch Time: " + str(las_file_date))

arcpy.env.extent = extent_fc

def pre_processing():
    """
    This function performs the preprocessing of the LiDAR las file required to create the CSV file. The output is a
    point feature class containing LiDAR points between the user defined starting height and the CHM offset, and a field
    for the Z value, CHM, and DTM.
    """

    print("\nPreprocessing\n")

    if os.path.exists(lidar_file_conversion):
        print("Deleting " + lidar_file_conversion)
        os.remove(lidar_file_conversion)

    if os.path.exists(lidar_clip_file):
        print("Deleting " + lidar_clip_file)
        os.remove(lidar_clip_file)

    las_file_clip = input_las_file.split(".las")[0] + "_clip" + ".las"
    if os.path.exists(las_file_clip):
        print("Deleting " + las_file_clip)
        os.remove(las_file_clip)

    print("Converting LAS file to an ArcGIS compatible LAS dataset file...")

    arcpy.conversion.ConvertLas(
        NEON_lidar_laz_file,
        lidar_folder,
        "SAME_AS_INPUT", '', "NO_COMPRESSION", None,
        lidar_file_conversion,
        "NO_FILES", None)

    print("Extracting LAS dataset file to study area...")
    # "_clip" is the suffix name added to the las file base name (used below in LASToMultipoint).
    arcpy.ddd.ExtractLas(
        lidar_file_conversion,
        lidar_folder,
        "",
        extent_fc, "PROCESS_EXTENT", '_clip', "MAINTAIN_VLR", "REARRANGE_POINTS",
        "COMPUTE_STATS",
        lidar_clip_file,
        "SAME_AS_INPUT")

    print("Converting LAS file to Multipoints...")
    # The extremely small point spacing parameter forces 1 point for each lidar point.
    arcpy.ddd.LASToMultipoint(
        las_file_clip,
        lidar_multipoint,
        1E-07, [], "ANY_RETURNS", None,
        output_proj,
        "las", 1, "NO_RECURSION")

    print("Adding Z value to Multipoints...")
    # Z_max, Z_min, Z_mean all return the same thing since it's a point.
    # Takes a long time. Maybe try adding z value to points after created below?
    arcpy.ddd.AddZInformation(lidar_multipoint, 'Z_max', 'NO_FILTER')

    print("Converting Multipoint to points...")  # Needed to extract by mask
    arcpy.management.FeatureToPoint(lidar_multipoint, lidar_multipoint_to_points, "CENTROID")

    print("Extracting DTM to Study Area...")
    # Only used for display purposes
    out_raster_dtm = arcpy.sa.ExtractByMask(NEON_DTM, extent_fc, "INSIDE")
    out_raster_dtm.save(output_dtm)

    print("Extracting CHM to Study Area...")
    # Only used for display purposes
    out_raster_chm = arcpy.sa.ExtractByMask(NEON_CHM, extent_fc, "INSIDE")
    out_raster_chm.save(output_chm)

    print("Extracting DTM to points...")
    # The extent comes up short of the points.
    arcpy.sa.ExtractValuesToPoints(lidar_multipoint_to_points, NEON_DTM, lidar_points_with_dtm_extraction, "", "VALUE_ONLY")
    arcpy.AlterField_management(lidar_points_with_dtm_extraction, "RASTERVALU", "dtm_extraction", "dtm_extraction")

    print("Extracting CHM to points...")
    # The extent comes up short of the points.
    arcpy.sa.ExtractValuesToPoints(lidar_points_with_dtm_extraction, NEON_CHM, lidar_points_with_z_dtm_and_chm, "", "VALUE_ONLY")
    arcpy.AlterField_management(lidar_points_with_z_dtm_and_chm, "RASTERVALU", "chm_extraction", "chm_extraction")

    print("\nCalculating height of each point from ground (DTM).")
    print("Points less than " + str(starting_height) + " meter from the ground will be dropped.\n")
    arcpy.AddField_management(lidar_points_with_z_dtm_and_chm, "height_from_ground", "LONG")

    with arcpy.da.UpdateCursor(lidar_points_with_z_dtm_and_chm, ["Z_max", "dtm_extraction", "chm_extraction", "height_from_ground"]) as uc:
        for row in uc:
            chm = row[2]  # Max height = CHM
            height_from_ground = row[0] - row[1]  # Height from ground = Z value from LiDAR - DTM

            # Delete LiDAR Point errors (points greater than the threshold above the CHM)
            if height_from_ground > chm + max_chm_offset:
                offset = height_from_ground - chm
                print("Deleting point with height = " + str(height_from_ground) + ". CHM is " + str(chm) + ". This point is " + str(offset) + " above the CHM. Max CHM offset is " + str(max_chm_offset))
                uc.deleteRow()
            elif starting_height and height_from_ground <= starting_height:
                # Too many print statements
                #print("Deleting point with height = " + str(height_from_ground) + ". Starting height is " + str(starting_height))
                uc.deleteRow()
            else:
                if height_from_ground < 0:
                    height_from_ground = 0
                row[3] = height_from_ground
                uc.updateRow(row)


def create_csv():
    """
    This function creates a CSV file containing LiDAR point counts aggregated to a regular volumetric grid (X,Y,Z) like
    a Rubik's cube. This is required in order to create a voxel layer.
    The lidar point returns are not initially structured like this; they are scattered in every x,y,z direction.
    The values in the z_max and z_min variables are used to determine the elevation "slices" used.
    For example a z_min of 1000 would mean that the first layer of voxels would be at an elevation of 1000m,
    and a z_max of 2000 would mean that the last f voxels would be at an elevation of 2000m.
    """

    print("\nCreating CSV\n")

    # Set the max elevation slice to be used in the "cube".
    all_z_values = [i[0] for i in arcpy.da.SearchCursor(lidar_points_with_z_dtm_and_chm, "Z_max")]
    z_max = int(round(max(all_z_values), 1))

    # Set the min elevation slice to be used in the "cube".
    # Use DTM for the min z value instead of min Lidar Z value because of a lidar point error that was ~500m below the surface.
    # This resulted in excessive iterations.
    all_dtm_values = [i[0] for i in arcpy.da.SearchCursor(lidar_points_with_z_dtm_and_chm, "dtm_extraction")]
    z_min = int(min(all_dtm_values))

    print("Min Global Elevation (from DTM): " + str(z_min))
    print("Max Global Elevation (from LiDAR Z values after errors removed): " + str(z_max))
    print("These values are used to determine the starting and ending elevation 'slices'.")

    z_start = z_min
    z_end = z_min + voxel_size

    # Create a fishnet and fishnet points for each z slice containing a count of the point returns.
    tmp_fishnet_points_to_merge = []

    count = 1
    #while z_end <= 1210:  # For testing
    while z_end <= z_max + 1:
        print("\nElevation Slice: " + str(count) + "/" + str(z_max - z_min))
        print("Starting Elevation: " + str(z_start))
        print("Ending Elevation: " + str(z_end))

        print("Creating Fishnet for this elevation slice...")

        tmp_fishnet = os.path.join(tmp_gdb, "fishnet_" + str(z_start) + "_" + str(z_end))
        desc = arcpy.Describe(output_dtm)
        arcpy.CreateFishnet_management(tmp_fishnet, str(desc.extent.lowerLeft),
                                       str(desc.extent.XMin) + " " + str(desc.extent.YMax), voxel_size, voxel_size, None, None,
                                       str(desc.extent.upperRight), "NO_LABELS", "#", "POLYGON")

        arcpy.AddField_management(tmp_fishnet, "Z", "DOUBLE")
        arcpy.CalculateField_management(tmp_fishnet, "Z", z_start)

        tmp_fishnet_with_point_count = os.path.join(tmp_gdb, "fishnet_with_point_count_" + str(z_start) + "_" + str(z_end))
        expression = "Z_max >= " + str(z_start) + " And Z_max < " + str(z_end)

        print("Joining LiDAR points to fishnet using CONTAINS method...")

        input_point_layer = arcpy.MakeFeatureLayer_management(lidar_points_with_z_dtm_and_chm)
        arcpy.management.SelectLayerByAttribute(input_point_layer, "NEW_SELECTION", expression, None)
        arcpy.analysis.SpatialJoin(tmp_fishnet,
                                   input_point_layer,
                                   tmp_fishnet_with_point_count, "JOIN_ONE_TO_ONE", "KEEP_ALL",'',
                                   "CONTAINS", None, '')

        arcpy.AlterField_management(tmp_fishnet_with_point_count, "Join_Count", "point_returns", "point_returns")

        print("Converting fishnet to points...")
        tmp_fishnet_points = os.path.join(tmp_gdb, "fishnet_with_point_count_point_conversion_" + str(z_start) + "_" + str(z_end))
        arcpy.FeatureToPoint_management(tmp_fishnet_with_point_count, tmp_fishnet_points)

        arcpy.management.DefineProjection(tmp_fishnet_points, output_proj)

        print("Adding an X field...")
        arcpy.AddField_management(tmp_fishnet_points, "X", "Double")

        print("Adding a Y field...")
        arcpy.AddField_management(tmp_fishnet_points, "Y", "Double")

        print("Calculating X and Y Coordinates...")
        arcpy.management.CalculateGeometryAttributes(
            tmp_fishnet_points,
            "X POINT_X;Y POINT_Y", '', '', None, "SAME_AS_INPUT")

        tmp_fishnet_points_to_merge.append(tmp_fishnet_points)

        z_start += 1
        z_end += 1
        count += 1

    print("\nMerging Fishnet Points...")
    arcpy.Merge_management(tmp_fishnet_points_to_merge, merged_fishnet_points)
    arcpy.AddField_management(merged_fishnet_points, "Time", "LONG")
    arcpy.CalculateField_management(merged_fishnet_points, "Time", las_file_date)
    #arcpy.AddField_management(merged_fishnet_points, "Time", "TEXT")
    #arcpy.CalculateField_management(merged_fishnet_points, "Time", "'" + las_file_date + "'")

    print("Exporting to CSV...")
    output_csv_dir = os.path.dirname(output_csv)
    output_csv_name = os.path.basename(output_csv)
    arcpy.TableToTable_conversion(merged_fishnet_points, output_csv_dir, output_csv_name)


def create_voxel():

    """ This function creates the NetCDF voxel layer from the CSV. The code is a modified version of the "Create a voxel
    layer from a CSV file" code available on the ESRI website:
    https://pro.arcgis.com/en/pro-app/latest/help/mapping/layer-properties/create-a-voxel-layer.htm """

    print("\nCreating Voxel\n")

    #Create a pandas dataframe and insert data from CSV/TEXT file
    dfPoints = pd.read_csv(output_csv)

    #Sort values to ensure they are in the correct order
    dfPoints = dfPoints.sort_values(by=['Time', 'Z','Y','X'])

    #Create domain for longitude/latitude
    #Each domain has unique values, no repeating numbers, and are sorted (to be monotonic)
    xDomain = np.sort(np.unique(dfPoints.iloc[:,10].values)) # 0th column contains x values
    print("X Domain Values:")
    print(xDomain)
    yDomain = np.sort(np.unique(dfPoints.iloc[:,11].values)) # 1st column contains y values
    print("Y Domain Values:")
    print(yDomain)
    zDomain = np.sort(np.unique(dfPoints.iloc[:,3].values)) # 2nd column contains z values
    print("Z Domain Values:")
    print(zDomain)
    tDomain = np.sort(np.unique(dfPoints.iloc[:,12].values)) # 3rd column contains t values
    print("T Domain Values:")
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
    #ncT = outDataSet.createVariable('t', np.float32, 't') # Create variable t
    ncT = outDataSet.createVariable('t', np.int32, 't') # Create variable t
    #ncT = outDataSet.createVariable('t', np.datetime64, 't') # Create variable t
    #ncT = outDataSet.createVariable('t', np.str_, 't') # Create variable t
    #ncT = outDataSet.createVariable('t', np.unicode_, 't') # Create variable t

    #Create variable data with dimensions (t,z,y,x). The fill value is set to -99999 which are values to be ignored by client.
    ncData = outDataSet.createVariable('data',np.float32,('t','z','y','x'),fill_value = -99999)

    #Assign values1
    ncX[:] = xDomain[:]
    ncY[:] = yDomain[:]
    ncZ[:] = zDomain[:]
    ncT[:] = tDomain[:]

    #test = np.arange(900)
    #should be 900 points
    #The dataframe 'Data' column must be reshaped to match the dimension shapes and placed into the ncData variable
    ncData[:,:,:,:] = np.reshape(
        dfPoints['point_returns'].values,
        #test,
        (tDomain.shape[0],zDomain.shape[0],yDomain.shape[0],xDomain.shape[0])
        )

    #Assign variable attributes

    ncData.long_name = "LiDAR Point Return Count"

    ncZ.positive = 'up'

    ncX.standard_name = 'projection_x_coordinate'
    ncX.units = 'm'

    ncY.standard_name = 'projection_y_coordinate'
    ncY.units = 'm'

    ncT.units = 'seconds since 1970-01-01 00:00:00'

    #Assign global attribute. This attribute is to assign a coordinate system.
    outDataSet.esri_pe_string = 'PROJCS["WGS_1984_UTM_Zone_11N",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",-117.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]],VERTCS["unknown",VDATUM["unknown"],PARAMETER["Vertical_Shift",0.0],PARAMETER["Direction",1.0],UNIT["Meter",1.0]]',

    outDataSet.close()

pre_processing()
create_csv()
create_voxel()

end_script = datetime.now()
print("\nEnd Time: " + str(end_script))
print("Duration: " + str(end_script - start_script))

