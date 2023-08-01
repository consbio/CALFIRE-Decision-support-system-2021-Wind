########################################################################################################################
# Author: Mike Gough
# Date created: 07/15/2023
# Python Version: 3.x
# Description:
# Merges the plot centroids and extracts CFO raster data for each centroid. Also intersects with the plots to get the
# plot id.
########################################################################################################################

import arcpy
arcpy.env.overwriteOutput = True

point_gdb = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\veg_structure\get_rows_and_colums\points.gdb"
plots = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\veg_structure\veg_structure.gdb\NEON_SOAP_veg_structure_2019_06_mbp_project"

raster_dirs = [
    r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\CFO_Data_Processing\Data\Outputs\CFO_Mosaics_Clipped_To_SSN_Summer_2020",
    r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\CFO_Data_Processing\Data\Outputs\CFO_Green_Vegetation_Clipped_To_SSN_Spring_2020",
    r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Text_File_to_GeoTiff\geotiff\b_lt_0"
]

merged_points = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\veg_structure\veg_structure.gdb\NEON_SOAP_merged_points_in_cells"
output_points = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\vegetation_structure\vegetation_structure.gdb\NEON_SOAP_veg_structure_plots_cell_centroids"

arcpy.env.workspace = point_gdb

points = arcpy.ListFeatureClasses()

arcpy.Merge_management(points, merged_points)

arcpy.Intersect_analysis([merged_points, plots], output_points)

for raster_dir in raster_dirs:
    arcpy.env.workspace = raster_dir
    rasters = arcpy.ListRasters()
    arcpy.sa.ExtractMultiValuesToPoints(output_points,rasters)

