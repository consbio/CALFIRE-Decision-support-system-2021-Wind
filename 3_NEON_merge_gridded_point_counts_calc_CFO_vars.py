########################################################################################################################
# File name: 3_NEON_merge_gridded_point_counts_calc_CFO_vars.py
# Author: Mike Gough
# Date created: 08/18/2023
# Python Version: 3.x
# Description:
# Combines the LiDAR point count data from a set of CSV files containing ALL returns (v8) with a version containing just
# the vegetation returns (v9) in order perform the "green" calculation below (comparable to CFO green vegetation, but
# using all veg returns, not just "green" vegetation since there is no way to get at that from the LiDAR data).
# This script also calculates canopy cover and ladder fuels using the equations below. All equations were provided
# by Dr. Joseph Werne with NorthWest Research Associates and are designed to match the calculations used by the
# California Forest Observatory (CFO) to make these same calculations:
#
#canopy cover = (total counts between canopy height and 5m above ground ) / (total count)
#ladder fuels = (total counts between 4m and 1m) / (total counts below 4m)
#green = (total counts in spectral range reflected from photosynthetic vegetation) / (total counts)

########################################################################################################################
import arcpy
import os
arcpy.env.overwriteOutput = True

version_with_all_points = "v8"
version_with_veg_points = "v9veg"  # Currently being generated
output_version = "v1"

input_gdb = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\Outputs.gdb"
intermediate_gdb = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Volume\Combine_Lidar_Point_Counts\Combine_Lidar_Point_Counts.gdb"
output_fc = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\CFO_Calcs_From_LiDAR\Merged_Outputs_With_CFO_Calculations.gdb\lidar_gridded_point_counts_merged_with_cfo_calcs_" + output_version
output_raster_canopy_cover = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\CFO_Calcs_From_LiDAR\GeoTiff" + os.sep + output_version + os.sep + "neon_soap_lidar_to_canopy_cover_" + output_version + ".tif"
output_raster_ladder_fuels = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\CFO_Calcs_From_LiDAR\GeoTiff" + os.sep + output_version + os.sep + "neon_soap_lidar_to_ladder_fuels_" + output_version + ".tif"
output_raster_green_veg = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\CFO_Calcs_From_LiDAR\GeoTiff" + os.sep + output_version + os.sep + "neon_soap_lidar_to_green_veg_" + output_version + ".tif"

all_points_merge_fc = os.path.join(intermediate_gdb, "lidar_gridded_point_counts_tile_fc_merge_" + version_with_all_points)
veg_points_merge_fc = os.path.join(intermediate_gdb, "lidar_gridded_point_counts_tile_fc_merge_" + version_with_veg_points)


def add_neon_tile_id():
    """ Adds a 'NEON tile ID' field to the individual tile-based outputs """

    print("Add a 'NEON tile ID' field and 'NEON tile ID and GRID ID' field to the individual tile-based outputs")

    arcpy.env.workspace = input_gdb
    all_points_fc_list = arcpy.ListFeatureClasses("*" + version_with_all_points)
    veg_points_fc_list = arcpy.ListFeatureClasses("*" + version_with_veg_points)

    for fc_list in [all_points_fc_list, veg_points_fc_list]:

        for tile_fc in fc_list:
            NEON_tile_id_parts = tile_fc.split("tile_")[-1].split("_")
            NEON_tile_id = str("_".join(NEON_tile_id_parts[:5]))
            arcpy.AddField_management(tile_fc, "NEON_tile_id", "TEXT")
            print(NEON_tile_id)
            arcpy.CalculateField_management(tile_fc, "NEON_tile_id", "'" + NEON_tile_id + "'")
            arcpy.AddField_management(tile_fc, "NEON_tile_and_grid_id", "TEXT")
            expression = "!NEON_tile_id!" + " + '_' + " + "str(!grid_id!)"
            arcpy.CalculateField_management(tile_fc, "NEON_tile_and_grid_id", expression)


def merge():
    """ Merges individual NEON tile feature classes into single FC for both ALL and just VEG versions."""

    print("Merge individual NEON tile feature classes into single FC for both ALL and just VEG versions.")
    arcpy.env.workspace = input_gdb

    # Merge output gridded point features classes that have ALL LiDAR points
    all_points_fc_list = arcpy.ListFeatureClasses("*" + version_with_all_points)

    #for all_point_fc in all_points_fc_list:
    #    NEON_tile_id_parts = all_point_fc.split("tile_")[-1].split("_")
    #    NEON_tile_id = "_".join(NEON_tile_id_parts[:5])
    #    print(NEON_tile_id)
    #    arcpy.AddField_management(all_point_fc, "NEON_tile_id", "TEXT")
    #    arcpy.CalculateField_management(all_point_fc, "NEON_tile_id", NEON_tile_id)

    arcpy.Merge_management(all_points_fc_list, all_points_merge_fc)

    # Merge output gridded point features classes that have just VEG LiDAR points
    veg_points_fc_list = arcpy.ListFeatureClasses("*" + version_with_veg_points)
    arcpy.Merge_management(veg_points_fc_list, veg_points_merge_fc)


def add_total_count_field():
    """ Creates a field containing a total count of all LiDAR point returns. """
    
    print("Add total count field if one is not already present")
    arcpy.env.workspace = intermediate_gdb
    merged_fc_list = arcpy.ListFeatureClasses()
    for merged_fc in merged_fc_list:
        # Just applied to v8, since a total count was not created for this version.
        fields = [field.name for field in arcpy.ListFields(merged_fc)]
        if "count_all" not in fields:
            # the veg version should have "count_all" in it, so skip this calculation.
            arcpy.AddField_management(merged_fc, "count_all", "LONG")
            fields_to_sum = ["count_min_1m", "count_1m_4m", "count_4m_5m", "count_5m_max"]
            fields_for_uc = fields_to_sum + ["count_all"]
            print(fields_for_uc)
            with arcpy.da.UpdateCursor(merged_fc, fields_for_uc) as uc:
                for row in uc:
                     row[4] = row[0] + row[1] + row[2] + row[3]
                     uc.updateRow(row)


def join_tables():
    """ Joins the VEG returns to the ALL returns feature class. """
    print("Joining veg returns and ALL returns feature classes...")
    arcpy.Copy_management(all_points_merge_fc, output_fc)
    arcpy.JoinField_management(output_fc,"NEON_tile_and_grid_id",veg_points_merge_fc,"NEON_tile_and_grid_id")


def calc_cfo_metrics():
    """ Calculate CFO metrics as defined by Joe Werne. Multiplied by 100 to match CFO values"""
    print("Adding CFO fields...")
    cfo_fields = ["canopy_cover", "ladder_fuels", "green_veg"]
    existing_fields = [field.name for field in arcpy.ListFields(output_fc)]
    for cfo_field in cfo_fields:
        if cfo_field not in existing_fields:
            print("Adding field: " + cfo_field)
            arcpy.AddField_management(output_fc, cfo_field, "DOUBLE")

    print("Calculating CFO Metrics...")
    fields_for_uc = cfo_fields + ["count_5m_max", "count_all", "count_1m_4m", "count_min_1m", "count_all_1"]
    print("Fields for Update Cursor: " + str(fields_for_uc))
    with arcpy.da.UpdateCursor(output_fc, fields_for_uc) as uc:
        for row in uc:
            count_5m_max = float(row[3])
            count_all = float(row[4])
            count_1m_4m = float(row[5])
            count_min_1m = float(row[6])
            count_all_veg = float(row[7])

            count_lt_4m = count_min_1m + count_1m_4m

            # Skip if no LiDAR Point Returns at all...
            if count_all:
                # Canopy Cover
                row[0] = (count_5m_max / count_all) * 100
                # Ladder Fuels
                if count_lt_4m:
                    row[1] = (count_1m_4m / count_lt_4m) * 100 # I confirmed that in places where count_lt_4m is 0, count_1m_4m is also 0.
                # Green Vegetation
                row[2] = (count_all_veg / count_all) * 100

            uc.updateRow(row)


def export_cfo_calcs_to_raster():
    """ Exports the final CFO calculations to GeoTiffs. """
    print("Exporting CFO Calculations to GeoTiffs...")
    arcpy.PolygonToRaster_conversion(output_fc, "canopy_cover", output_raster_canopy_cover, "CELL_CENTER", "NONE", 10)
    arcpy.PolygonToRaster_conversion(output_fc, "ladder_fuels", output_raster_ladder_fuels, "CELL_CENTER", "NONE", 10)
    arcpy.PolygonToRaster_conversion(output_fc, "green_veg", output_raster_green_veg, "CELL_CENTER", "NONE", 10)


add_neon_tile_id()
merge()
add_total_count_field()
join_tables()
calc_cfo_metrics()
export_cfo_calcs_to_raster()
