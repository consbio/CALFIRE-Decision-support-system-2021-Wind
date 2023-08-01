########################################################################################################################
# File name: get_rows_and_cols.py
# Author: Mike Gough
# Date created: 07/15/2023
# Python Version: 3.x
# Description:
# Gets the row and column range for a sub-area within a larger raster dataset.
# For example, given the CFO raster, find the row and column range for the NEON SOAP study site.
# This was requested for Joe Werne to allow him to process the CFO data within the NEON SOAP study site.
########################################################################################################################
import math
import arcpy
import pandas as pd
import os
from osgeo import gdal
arcpy.env.overwriteOutput = True

input_raster_dir = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\veg_structure\geotiff_clips\v1"
output_csv = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\CSV\row_and_column_indices\row_and_column_indices_raster_cells_in_NEON_plots.csv"
point_gdb = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\veg_structure\get_rows_and_colums\points.gdb"
raster_gdb = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\veg_structure\get_rows_and_colums\rasters.gdb"

#raster_columns = 26577
#raster_rows = 35595

# CFO (SSN Study Area) raster extent
raster_resolution = 10

raster_left = 688080
raster_right = 953850
raster_top = 4275130
raster_bottom = 3919180

arcpy.env.workspace = input_raster_dir

# Only need to get the col and row values from one raster at each plot. Chose "gam005" and "Ns" arbitrarily.
rasters = arcpy.ListRasters("gam005*_Ns_*")
range_dict = {}

for raster in rasters:
    print(raster)
    basename = raster.split(".")[0]

    # Convert clipped raster to points and then back to a raster to get rid of NoData pixels
    points = os.path.join(point_gdb, basename)
    arcpy.conversion.RasterToPoint(raster, points, raster_field="Value")

    points_to_raster = os.path.join(raster_gdb, basename)
    with arcpy.EnvManager(snapRaster=raster):
        arcpy.conversion.PointToRaster(
            in_features=points,
            value_field="grid_code",
            out_rasterdataset=points_to_raster,
            cell_assignment="MOST_FREQUENT",
            priority_field="NONE",
            cellsize=raster_resolution,
            build_rat="BUILD"
        )

    extent = arcpy.Describe(points_to_raster).extent
    boundary_left = extent.XMin
    boundary_right = extent.XMax
    boundary_top = extent.YMax
    boundary_bottom = extent.YMin

    column_start = int(math.floor((boundary_left - raster_left)/raster_resolution))
    column_end = int(math.ceil((boundary_right - raster_left)/raster_resolution)) - 1

    row_start = int(math.floor((raster_top - boundary_top)/raster_resolution))
    row_end = int(math.ceil((raster_top - boundary_bottom)/raster_resolution)) - 1

    column_range = str(column_start) + " - " + str(column_end)
    print("Column Range: " + column_range)

    row_range = str(row_start) + " - " + str(row_end)
    print("Row Range: " + row_range)

    #print("Column Start: " + str(column_start))

    xllcorner = raster_left + (column_start * 10)
    yllcorner = raster_top - (row_end * 10)

    #print("xllcorner: " + str(xllcorner))
    #print("yllcorner: " + str(yllcorner))

    # Pattern: gam005_0_fG_clip_SOAP_001.tif
    NEON_plot = raster.split("clip_")[-1].replace(".tif", "")
    range_dict[NEON_plot] = {}
    range_dict[NEON_plot]["cfo_ssn_row_range"]=row_range
    range_dict[NEON_plot]["cfo_ssn_column_range"]=column_range

print("\nWrite to CSV...")
df = pd.DataFrame.from_dict(range_dict, orient="index")
df.index.name = "plotID"
df.to_csv(output_csv, index=True, header=True)


