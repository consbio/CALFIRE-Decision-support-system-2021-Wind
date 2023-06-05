########################################################################################################################
# File name: NEON_count_LiDAR_point_returns.py
# Author: Mike Gough
# Date created: 06/01/2023
# Python Version: 2.7
# Description:
# Counts the number of LiDAR point returns for a set of regular vertical height intervals measured from the ground(DTM).
# The output is a fishnet (vector based grid) with a field containing a point count for each interval
# For example, the field count_0m_1m stores the number of point returns recorded between 0m and 1m from the ground for
# each polygon grid.
#######################################################################################################################
import arcpy
import os

arcpy.env.overwriteOutput = True

extent_fc = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Inputs\Extents\Extents.gdb\Extent_1_Tower_Location"
extent_name = extent_fc.split(os.sep)[-1]

fishnet_input_points_extent = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Inputs\Volume\Vector_Fishnets.gdb\Fishnet_LiDAR_Point_" + extent_name

#input_points_with_z = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermedate\Volume\Volume.gdb\NEON_D17_SOAP_Las_to_Multipoint_Subset_Small_Subset_to_Point"
input_points_with_z_and_height_from_ground = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Volume\Volume.gdb\Lidar_Points_with_Elevation_" + extent_name

dtm = r"\\loxodonta\gis\Source_Data\environment\region\NEON_SITES\SOAP\Elevation_LiDAR\2021\07\NEON_lidar-elev\NEON.D17.SOAP.DP3.30024.001.2021-07.basic.20230601T181117Z.RELEASE-2023\NEON_D17_SOAP_DP3_298000_4100000_DTM.tif"
tmp_gdb = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Scratch\Scratch.gdb"
intermediate_gdb = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Volume\Volume.gdb"
intermediate_folder = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Volume"
output_fc = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\Outputs.gdb\segmented_heights_above_ground_" + extent_name

height_interval = 1

# arcpy.management.CreateFishnet(output_grid, "297023.121 4105071.399", "297023.121 4105081.399", 1, 1, None, None, "297218.474 4105255.21", "NO_LABELS", '297023.121 4105071.399 297218.474 4105255.21 PROJCS["WGS_1984_UTM_Zone_11N",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",-117.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]', "POLYGON")
arcpy.env.extent = extent_fc


def pre_processing():
    NEON_lidar_laz_file = r"\\loxodonta\gis\Source_Data\environment\region\NEON_SITES\SOAP\Discrete_return_LiDAR_point_cloud\2019\06\NEON_lidar-point-cloud-line\NEON.D17.SOAP.DP1.30003.001.2019-06.basic.20230523T232633Z.RELEASE-2023\NEON_D17_SOAP_DP1_298000_4100000_classified_point_cloud_colorized.laz"
    NEON_DTM = r"\\loxodonta\gis\Source_Data\environment\region\NEON_SITES\SOAP\Elevation_LiDAR\2021\07\NEON_lidar-elev\NEON.D17.SOAP.DP3.30024.001.2021-07.basic.20230601T181117Z.RELEASE-2023\NEON_D17_SOAP_DP3_298000_4100000_DTM.tif"

    print("Converting LAS file...")

    lidar_folder = intermediate_folder + os.sep + "Lidar"
    lidar_file_conversion_name = NEON_lidar_laz_file.split(os.sep)[-1].split(".")[0] + "_LiDAR_LAS_File_Conversion.lasd"
    lidar_file_conversion = lidar_folder + os.sep + lidar_file_conversion_name

    arcpy.conversion.ConvertLas(
        NEON_lidar_laz_file,
        lidar_folder,
        "SAME_AS_INPUT", '', "NO_COMPRESSION", None,
        lidar_file_conversion,
        "NO_FILES", None)

    lidar_clip_file = lidar_folder + os.sep + lidar_file_conversion_name.split(".")[0] + "_Clip" + ".lasd"
    print("Extracting LAS file to study area...")
    arcpy.ddd.ExtractLas(
        lidar_file_conversion,
        lidar_folder,
        "",
        extent_fc, "PROCESS_EXTENT", '', "MAINTAIN_VLR", "REARRANGE_POINTS",
        "COMPUTE_STATS",
        lidar_clip_file,
        "SAME_AS_INPUT")

    print("Converting LAS file to Multipoints...")
    input_las_file = lidar_folder + os.sep + NEON_lidar_laz_file.split(os.sep)[-1].replace(".laz", ".las")
    lidar_multipoint = intermediate_gdb + os.sep + "lidar_multipoint"
    arcpy.ddd.LASToMultipoint(
        input_las_file,
        lidar_multipoint,
        1E-07, [], "ANY_RETURNS", None,
        'PROJCS["WGS_1984_UTM_Zone_11N",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",-117.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]],VERTCS["unknown",VDATUM["unknown"],PARAMETER["Vertical_Shift",0.0],PARAMETER["Direction",1.0],UNIT["Meter",1.0]]',
        "las", 1, "NO_RECURSION")

    print("Adding Z value to Multipoints...")
    arcpy.ddd.AddZInformation(lidar_multipoint, 'Z_max', 'NO_FILTER')

    print("Converting LAS file to points...")
    lidar_multipoint_to_points = intermediate_gdb + os.sep + "lidar_multipoint_to_points"
    arcpy.management.FeatureToPoint(lidar_multipoint, lidar_multipoint_to_points, "CENTROID")

    print("Extracting DTM to Study Area...")
    output_dtm = intermediate_gdb + os.sep + "NEON_DTM_Clipped_" + extent_fc.split(os.sep)[-1]
    out_raster = arcpy.sa.ExtractByMask(NEON_DTM, extent_fc, "INSIDE")
    out_raster.save(output_dtm)

    print("Extracting DTM to points...")
    # The extent comes up short of the points.
    arcpy.sa.ExtractValuesToPoints(lidar_multipoint_to_points, NEON_DTM, input_points_with_z_and_height_from_ground, "", "VALUE_ONLY")

    print("Creating Fishnet...")
    desc = arcpy.Describe(output_dtm)
    arcpy.CreateFishnet_management(fishnet_input_points_extent, str(desc.extent.lowerLeft),
                                   str(desc.extent.XMin) + " " + str(desc.extent.YMax), "1", "1", None, None,
                                   str(desc.extent.upperRight), "NO_LABELS", "#", "POLYGON")


def count_point_returns():

    arcpy.CopyFeatures_management(fishnet_input_points_extent, output_fc)

    #input_points_with_z_and_height_from_ground = intermediate_gdb + os.sep + "NEON_D17_SOAP_Las_to_Multipoint_Subset_Small_Subset_to_Point_Extract_DTM"

    #arcpy.sa.ExtractValuesToPoints(input_points_with_z, dtm, input_points_with_z_and_height_from_ground, "NONE", "VALUE_ONLY")

    arcpy.AddField_management(input_points_with_z_and_height_from_ground, "height_from_ground", "DOUBLE")

    print("Calculating height of each point from ground (DTM)")
    with arcpy.da.UpdateCursor(input_points_with_z_and_height_from_ground, ["Z_max", "RASTERVALU", "height_from_ground"]) as uc:
        for row in uc:
            height_from_ground = row[0] - row[1]
            if height_from_ground < 0:
                height_from_ground = 0
            row[2] = height_from_ground
            uc.updateRow(row)

    all_h_values = [i[0] for i in arcpy.da.SearchCursor(input_points_with_z_and_height_from_ground,"height_from_ground")]

    h_min = int(min(all_h_values))
    h_max = int(round(max(all_h_values), 1))

    print("Min Height: " + str(h_min))
    print("Max Height: " + str(h_max))

    h_start = h_min
    h_end = h_min + height_interval

    while h_end <= h_max + 1:

        print("Starting Height: " + str(h_start))
        print("Ending Height: " + str(h_end))

        output_spatial_join_fc = tmp_gdb + os.sep + "spatial_join_" + str(h_start) + "_" + str(h_end)
        expression = "height_from_ground >= " + str(h_start) + " And height_from_ground < " + str(h_end)
        #print(expression)
        input_point_layer = arcpy.MakeFeatureLayer_management(input_points_with_z_and_height_from_ground)
        arcpy.management.SelectLayerByAttribute(input_point_layer, "NEW_SELECTION", expression, None)

        print("Performing spatial join to get a count of the number LiDAR point returns within this height segment for each polygon...")
        arcpy.analysis.SpatialJoin(output_fc, input_point_layer, output_spatial_join_fc, "JOIN_ONE_TO_ONE", "KEEP_ALL", '#' "CONTAINS", None, '')

        print("Creating a dictionary to store these returns (OBJECTID[point_count])")
        count_dict = {}

        with arcpy.da.SearchCursor(output_spatial_join_fc, ["TARGET_FID", "Join_Count"]) as sc:
            for row in sc:
                count_dict[row[0]] = row[1]

        print("Adding point count for this height segment to the output Feature Class...\n")
        field_name = "count_" + str(h_start) + "m_" + str(h_end) + "m"
        arcpy.AddField_management(output_fc, field_name, "LONG")
        with arcpy.da.UpdateCursor(output_fc, ["OBJECTID", field_name]) as uc:
            for row in uc:
                row[1] = count_dict[row[0]]
                uc.updateRow(row)

        h_start += 1
        h_end += 1


pre_processing()
count_point_returns()
