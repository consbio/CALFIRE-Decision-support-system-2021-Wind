########################################################################################################################
# File name: 1_NEON_lidar_to_gridded_point_count_at_h.py
# Author: Mike Gough
# Date created: 06/01/2023
# Python Version: 3.x
# Description:
# Counts the number of LiDAR point returns for a set of regular vertical height intervals measured from the ground(DTM).
# The output is a fishnet (vector based grid) with a field containing a point count for each interval
# For example, the field count_0m_1m stores the number of point returns recorded between 0m and 1m from the ground for
# each polygon grid.
########################################################################################################################
import arcpy
import os
from datetime import datetime

arcpy.env.overwriteOutput = True

version = "v3_test_ground_returns"
height_interval = 1
output_resolution = 1
max_chm_offset = 10  # Max Canopy Height Offset allowed. LiDAR point returns greater than this distance above the CHM will be considered errors and will be deleted (units are m).


extent_fc = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Inputs\Extents\Extents.gdb\extent_1_tower_location"
#extent_fc = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Inputs\Extents\Extents.gdb\extent_1_tower_location_utm_zone_10n"
extent_name = extent_fc.split(os.sep)[-1]


NEON_lidar_laz_file = r"\\loxodonta\gis\Source_Data\environment\region\NEON_SITES\SOAP\Discrete_return_LiDAR_point_cloud\2019\06\NEON_lidar-point-cloud-line\NEON.D17.SOAP.DP1.30003.001.2019-06.basic.20230523T232633Z.RELEASE-2023\NEON_D17_SOAP_DP1_298000_4100000_classified_point_cloud_colorized.laz"
NEON_DTM = r"\\loxodonta\gis\Source_Data\environment\region\NEON_SITES\SOAP\Elevation_LiDAR\2021\07\NEON_lidar-elev\NEON.D17.SOAP.DP3.30024.001.2021-07.basic.20230601T181117Z.RELEASE-2023\NEON_D17_SOAP_DP3_298000_4100000_DTM.tif"
NEON_CHM = r"\\loxodonta\gis\Source_Data\environment\region\NEON_SITES\SOAP\Ecosystem_structure\2019\06\NEON_struct-ecosystem\NEON.D17.SOAP.DP3.30015.001.2019-06.basic.20230524T172838Z.RELEASE-2023\NEON_D17_SOAP_DP3_298000_4100000_CHM.tif"
snap_grid = r"\\loxodonta\gis\Projects\CALFIRE_Decision_support_system_2021\Workspaces\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\CFO_Data_Processing\Data\Outputs\CFO_Green_Vegetation_By_County_Spring_2020\Fresno-County-California-Vegetation-GreenVegetation-2020-Spring-00010m.tif"

tmp_gdb = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Scratch\Scratch.gdb"
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
output_snap_grid = os.path.join(intermediate_gdb, "SNAP_GRID_clipped_" + extent_name + "_" + version)
lidar_points_with_dtm_extraction = os.path.join(tmp_gdb, "lidar_points_with_dtm_extraction_" + extent_name + "_" + version)
lidar_points_with_z_dtm_and_chm = os.path.join(intermediate_gdb, "lidar_points_with_z_dtm_and_chm_" + extent_name + "_" + version)

lidar_points_with_z_dtm_and_chm = os.path.join(intermediate_gdb, "lidar_points_with_elevation_" + extent_name + "_" + version)
fishnet_input_points_extent = os.path.join(intermediate_gdb, "fishnet_lidar_point_" + extent_name + "_" + version)

# Output Feature Class

output_fc = "_".join([r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\Outputs.gdb\segmented_heights_above_ground", extent_name, version])

start_script = datetime.now()
print("\nStart Time: " + str(start_script))

output_proj = arcpy.Describe(extent_fc).spatialReference

arcpy.env.extent = extent_fc
arcpy.env.outputCoordinateSystem = output_proj

print(fishnet_input_points_extent)
print(output_fc)

def pre_processing():

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

    print("Converting LAS file...")

    arcpy.conversion.ConvertLas(
        NEON_lidar_laz_file,
        lidar_folder,
        "SAME_AS_INPUT", '', "NO_COMPRESSION", None,
        lidar_file_conversion,
        "NO_FILES", None)

    print("Extracting LAS file to study area (and projecting LAS)...")
    lidar_clip_file = lidar_folder + os.sep + lidar_file_conversion_name.split(".")[0] + "_Clip" + ".lasd"
    with arcpy.EnvManager(outputCoordinateSystem=output_proj):
        arcpy.ddd.ExtractLas(
            lidar_file_conversion,
            lidar_folder,
            "",
            extent_fc, "PROCESS_EXTENT", '', "MAINTAIN_VLR", "REARRANGE_POINTS",
            "COMPUTE_STATS",
            lidar_clip_file,
            "SAME_AS_INPUT")

    print("Converting LAS file to Multipoints...")
    arcpy.ddd.LASToMultipoint(
        input_las_file,
        lidar_multipoint,
        1E-07, [], "ANY_RETURNS", None,
        #'PROJCS["WGS_1984_UTM_Zone_11N",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",-117.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]],VERTCS["unknown",VDATUM["unknown"],PARAMETER["Vertical_Shift",0.0],PARAMETER["Direction",1.0],UNIT["Meter",1.0]]',
        #output_proj,
        "",
        "las", 1, "NO_RECURSION")

    print("Adding Z value to Multipoints...")
    arcpy.ddd.AddZInformation(lidar_multipoint, 'Z_max', 'NO_FILTER')

    print("Converting Multipoint to points...")  # Needed to extract by mask
    arcpy.management.FeatureToPoint(lidar_multipoint, lidar_multipoint_to_points, "CENTROID")

    print("Extracting SNAP Grid to Study Area...")
    out_raster_snap_grid = arcpy.sa.ExtractByMask(snap_grid, extent_fc, "INSIDE")
    out_raster_snap_grid.save(output_snap_grid)

    print("Extracting DTM to Study Area...")
    out_raster_dtm = arcpy.sa.ExtractByMask(NEON_DTM, extent_fc, "INSIDE")
    out_raster_dtm.save(output_dtm)

    print("Extracting CHM to Study Area...")
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

    print("Creating Fishnet...")
    desc = arcpy.Describe(output_snap_grid)
    arcpy.CreateFishnet_management(fishnet_input_points_extent, str(desc.extent.lowerLeft),
                                   str(desc.extent.XMin) + " " + str(desc.extent.YMax), output_resolution, output_resolution, None, None,
                                   str(desc.extent.upperRight), "NO_LABELS", "#", "POLYGON")

    arcpy.management.DefineProjection(fishnet_input_points_extent, output_proj)

def count_point_returns():

    arcpy.CopyFeatures_management(fishnet_input_points_extent, output_fc)

    #lidar_points_with_z_dtm_and_chm = intermediate_gdb + os.sep + "NEON_D17_SOAP_Las_to_Multipoint_Subset_Small_Subset_to_Point_Extract_DTM"
    #arcpy.sa.ExtractValuesToPoints(input_points_with_z, dtm, lidar_points_with_z_dtm_and_chm, "NONE", "VALUE_ONLY")

    arcpy.AddField_management(lidar_points_with_z_dtm_and_chm, "height_from_ground", "DOUBLE")

    print("Calculating height of each point from ground (DTM)")
    with arcpy.da.UpdateCursor(lidar_points_with_z_dtm_and_chm, ["Z_max", "dtm_extraction", "chm_extraction", "height_from_ground",]) as uc:
        for row in uc:
            max_possible_height = row[2]  # Max height = CHM
            height_from_ground = row[0] - row[1]  # Height from ground = Z value from LiDAR - DTM

            # Delete LiDAR Point errors (points taller than the CHM)
            # This is causing problems. Ground points are being deleted because they are above the CHM which is 0.
            if height_from_ground > max_possible_height + max_chm_offset:
                print("Deleting point with height = " + str(height_from_ground) + ". Max height from CHM is " + str(max_possible_height))
                uc.deleteRow()
            else:
                if height_from_ground < 0:
                    height_from_ground = 0
                row[3] = height_from_ground
                uc.updateRow(row)

    all_h_values = [i[0] for i in arcpy.da.SearchCursor(lidar_points_with_z_dtm_and_chm, "height_from_ground")]

    h_min = int(min(all_h_values))
    h_max = int(round(max(all_h_values), 1))

    print("Min Height: " + str(h_min))
    print("Max Height: " + str(h_max))

    h_start = h_min
    h_end = h_min + height_interval

    while h_end <= h_max + 1:

        print("Starting Height: " + str(h_start))
        print("Ending Height: " + str(h_end))

        output_spatial_join_fc = tmp_gdb + os.sep + "spatial_join_" + str(h_start).replace('-', "neg") + "_" + str(h_end).replace("-", "neg")
        expression = "height_from_ground >= " + str(h_start) + " And height_from_ground < " + str(h_end)
        #print(expression)
        input_point_layer = arcpy.MakeFeatureLayer_management(lidar_points_with_z_dtm_and_chm)
        arcpy.management.SelectLayerByAttribute(input_point_layer, "NEW_SELECTION", expression, None)

        print("Performing spatial join to get a count of the number LiDAR point returns within this height segment for each polygon...")
        arcpy.analysis.SpatialJoin(output_fc, input_point_layer, output_spatial_join_fc, "JOIN_ONE_TO_ONE", "KEEP_ALL", '#' "CONTAINS", None, '')

        print("Creating a dictionary to store these returns (OBJECTID[point_count])")
        count_dict = {}

        with arcpy.da.SearchCursor(output_spatial_join_fc, ["TARGET_FID", "Join_Count"]) as sc:
            for row in sc:
                count_dict[row[0]] = row[1]

        print("Adding point count for this height segment to the output Feature Class...\n")
        field_name = "count_" + str(h_start).replace("-", "neg") + "m_" + str(h_end).replace("-", "neg") + "m"
        arcpy.AddField_management(output_fc, field_name, "LONG")
        with arcpy.da.UpdateCursor(output_fc, ["OBJECTID", field_name]) as uc:
            for row in uc:
                row[1] = count_dict[row[0]]
                uc.updateRow(row)

        h_start += 1
        h_end += 1


def post_processing():
    print("Adding GRID ID to final output...")
    tmp_points = os.path.join(tmp_gdb, "tmp_points")
    tmp_points_with_dtm = os.path.join(tmp_gdb, "tmp_points_with_dtm")
    arcpy.AddField_management(output_fc, "GRID_ID", "LONG")
    arcpy.CalculateField_management(output_fc, "GRID_ID", "!OBJECTID!")

    print("Adding DTM to final output...")
    arcpy.FeatureToPoint_management(output_fc, tmp_points)
    arcpy.sa.ExtractValuesToPoints(tmp_points, NEON_DTM, tmp_points_with_dtm)
    arcpy.JoinField_management(output_fc,"GRID_ID",tmp_points_with_dtm, "GRID_ID", ["RASTERVALU"])
    arcpy.AlterField_management(output_fc, "RASTERVALU", "DTM_Extraction", "DTM_Extraction")

pre_processing()
count_point_returns()
post_processing()

end_script = datetime.now()
print("\nEnd Time: " + str(end_script))
print("Duration: " + str(end_script - start_script))
