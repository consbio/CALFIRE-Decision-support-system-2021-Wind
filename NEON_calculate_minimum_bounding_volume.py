import arcpy
import math
arcpy.env.overwriteOutput = True

input_fishnet_with_xy = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Inputs\Volume\Vector_Fishnets.gdb\Fishnet_LiDAR_Point_Extent_1_Tower_Location_For_Volume_Calculation_Join_XY"
input_point_with_z = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Volume\Volume.gdb\Lidar_Points_with_Elevation_Extent_1_Tower_Location"
intersect_fc = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Volume\Volume.gdb\LiDAR_Points_with_X_Y_Z_Index"
output_fc = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\Outputs.gdb\minimum_bounding_volume_convex_hull"

height_interval = 1

print("Intersecting Fishnet and Input Points")
#arcpy.Intersect_analysis([input_point_with_z, input_fishnet_with_xy], intersect_fc, "ALL")
#arcpy.AddField_management(intersect_fc, "Z_Index")

if False:
    print("Adding integer based Z-Index")
    with arcpy.da.UpdateCursor(intersect_fc, ["Z_Max", "height_from_ground", "Z_Index"]) as uc:
        for row in uc:
            z_index = math.floor(row[1])
            if z_index == 0:
                uc.deleteRow()
            else:
                row[2] = z_index
                uc.updateRow(row)


arcpy.ddd.MinimumBoundingVolume(intersect_fc, "Z_Max", output_fc, "CONVEX_HULL", "LIST", "POINT_X;POINT_Y;Z_Index", "MBV_FIELDS")


