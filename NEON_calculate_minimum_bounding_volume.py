########################################################################################################################
# File name: NEON_calculate_minimum_volume_geometery.py
# Author: Mike Gough
# Date created: 06/05/2023
# Python Version: 2.7
# Description:
# Calculates the minimum bounding volume for a set of pre-processed LiDAR point returns.
# Requires a fishnet with X,Y coordinates of the centroid (spatial join with fishnet labels), and a set of input points
# with a z value and height_above_ground.
########################################################################################################################

import arcpy
import math
arcpy.env.overwriteOutput = True

input_fishnet_with_xy = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Inputs\Volume\Vector_Fishnets.gdb\Fishnet_LiDAR_Point_Extent_1_Tower_Location_For_Volume_Calculation_Join_XY"
input_point_with_z = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Volume\Volume.gdb\Lidar_Points_with_Elevation_Extent_1_Tower_Location"
intersect_fc = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Volume\Volume.gdb\LiDAR_Points_with_X_Y_Z_Index"

#output_fc = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\Outputs.gdb\minimum_bounding_volume_envelope"
#output_fc = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\Outputs.gdb\minimum_bounding_volume_sphere"
#output_fc = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\Outputs.gdb\minimum_bounding_volume_convex_hull"
output_fc = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\Outputs.gdb\minimum_bounding_volume_concave_hull"

print("Intersecting Fishnet and Input Points")
arcpy.Intersect_analysis([input_point_with_z, input_fishnet_with_xy], intersect_fc, "ALL")
arcpy.AddField_management(intersect_fc, "Z_Index")

print("Removing ground points. Adding integer based Z-Index")
with arcpy.da.UpdateCursor(intersect_fc, ["Z_Max", "height_from_ground", "Z_Index"]) as uc:
    for row in uc:
        z_index = math.floor(row[1])
        if z_index == 0:
            uc.deleteRow()
        else:
            row[2] = z_index
            uc.updateRow(row)

# Notes: Sphere creates volumes that extend beyond the 3D cubes.
arcpy.ddd.MinimumBoundingVolume(intersect_fc, "Z_Max", output_fc, "CONCAVE_HULL", "LIST", "POINT_X;POINT_Y;Z_Index", "MBV_FIELDS")

arcpy.AddField_management(output_fc, "MBV_Percent", "DOUBLE")

with arcpy.da.UpdateCursor(output_fc, ["MBV_Volume", "MBV_Percent"]) as uc:
    for row in uc:
        volumetric_percent = row[0] / 1 * 100
        row[1] = volumetric_percent
        uc.updateRow(row)



