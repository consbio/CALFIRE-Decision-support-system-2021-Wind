# Author: Mike Gough
# Date created: 07/26/2023
# Python Version: 3.x
# Description:
# This script compares the calculated tree and shrub density results produced by Joe Werne with the in-situ
# counts in the vegetation sampling data produced by NEON: https://data.neonscience.org/data-products/DP1.10098.001
# The output is a CSV file that contains the differences between these values for each value of gamma.
########################################################################################################################

import glob
import os
import csv
import arcpy
import numpy as np
import pandas as pd
arcpy.env.overwriteOutput = True

# Input vars & dirs

#version = "v1"
#veg_structure_dir = r"\\loxodonta\gis\Source_Data\environment\region\NEON_SITES\SOAP\Vegetation_structure\2019\06\NEON_struct-plant\NEON.D17.SOAP.DP1.10098.001.2019-06.basic.20230127T120753Z.RELEASE-2023"

#version = "v2" # Note: not a lot of data. Only 1 plot and 41 apparent individuals. Did not use the data from 03/2020
#veg_structure_dir = r"\\loxodonta\gis\Source_Data\environment\region\NEON_SITES\SOAP\Vegetation_structure\2020\03\NEON_struct-plant\NEON.D17.SOAP.DP1.10098.001.2020-03.basic.20230127T120753Z.RELEASE-2023"

version = "v42" # Note: not a lot of data. Only 1 plot and 41 apparent individuals. Did not use v2.
veg_structure_dir = r"\\loxodonta\gis\Source_Data\environment\region\NEON_SITES\SOAP\Vegetation_structure\2019\06\NEON_struct-plant\NEON.D17.SOAP.DP1.10098.001.2019-06.basic.20230127T120753Z.RELEASE-2023"
input_geotiffs_dir = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Text_File_to_GeoTiff\v42\geotiff"

plot_filters = None
raster_filters = None
#plot_filters = ["SOAP_010"]
#raster_filters = ["gam010_0_Ns"]

# Intermediate dirs
intermediate_gdb_dir = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\veg_structure"
intermediate_gdb = os.path.join(intermediate_gdb_dir, version + ".gdb")
geotiffs_clip_dir = os.path.join(r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\veg_structure\geotiff_clips", version)
tmp_gdb = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\veg_structure\tmp\tmp.gdb"
tmp_dir = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\veg_structure\tmp"

if not arcpy.Exists(intermediate_gdb):
    arcpy.CreateFileGDB_management(intermediate_gdb_dir, version)

if not arcpy.Exists(geotiffs_clip_dir):
    os.mkdir(geotiffs_clip_dir)

NEON_data_date = veg_structure_dir.split("NEON.")[-1].split(".")[5].replace("-", "_")

# Output file:
output_csv = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\CSV\veg_structure_validation\veg_structure_validation_" + version + "_" + NEON_data_date + ".csv"

# CRS Info
arcpy.env.workspace = input_geotiffs_dir
veg_structure_crs = arcpy.SpatialReference(text='PROJCS["WGS_1984_UTM_Zone_11N",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",-117.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]];-5120900 -9998100 450445547.391054;-100000 10000;-100000 10000;0.001;0.001;0.001;IsHighPrecision')
raster_crs = arcpy.Describe(arcpy.ListRasters()[0]).spatialReference

# Derived paths
base_name = "NEON_SOAP_veg_structure_" + NEON_data_date

plot_csv = glob.glob(veg_structure_dir + os.sep + "*perplotperyear*.csv")[0]
apparent_individual_csv = glob.glob(veg_structure_dir + os.sep + "*apparentindividual*.csv")[0]

points = os.path.join(intermediate_gdb, base_name + "_sampling_sites")
buffer = os.path.join(intermediate_gdb, base_name + "_buffer")
mbp = os.path.join(intermediate_gdb, base_name + "_mbp")
mbp_project = os.path.join(intermediate_gdb, base_name + "_mbp_project")

apparent_individual_table_name = base_name + "_apparent_individual"
apparent_individual_table = os.path.join(intermediate_gdb, apparent_individual_table_name)

shrubs_table_name = base_name + "_shrub"
shrubs_table = os.path.join(intermediate_gdb, shrubs_table_name)

trees_table_name = base_name + "_trees"
trees_table = os.path.join(intermediate_gdb, trees_table_name)

summary_statistics_shrubs_i = os.path.join(intermediate_gdb, base_name + "_shrubs_summary_stats_individual")
summary_statistics_trees_i = os.path.join(intermediate_gdb, base_name + "_trees_summary_stats_individual")

summary_statistics_shrubs = os.path.join(intermediate_gdb, base_name + "_shrubs_summary_stats")
summary_statistics_trees = os.path.join(intermediate_gdb, base_name + "_trees_summary_stats")


plot_dict = {}
with open(plot_csv, "r") as csvfile:
    reader = csv.reader(csvfile)
    next(reader)
    for row in reader:
        plot_id = row[6]
        if (plot_filters and plot_id in plot_filters) or (not plot_filters):
            x = row[13]
            y = row[14]
            plot_dict[plot_id] = {}
            plot_dict[plot_id]["x"] = x
            plot_dict[plot_id]["y"] = y
            plot_dict[plot_id]["shrub_count"] = 0
            plot_dict[plot_id]["tree_count"] = 0

def create_plot_fc():

    print("\nCreate NEON sampling plot polygon...")

    arcpy.management.XYTableToPoint(
        in_table=plot_csv,
        out_feature_class=points,
        x_field="easting",
        y_field="northing",
        z_field=None,
        coordinate_system=veg_structure_crs
    )

    arcpy.analysis.Buffer(
        in_features=points,
        out_feature_class=buffer,
        buffer_distance_or_field="10 Meters",
        line_side="FULL",
        line_end_type="ROUND",
        dissolve_option="NONE",
        dissolve_field=None,
        method="PLANAR"
    )

    arcpy.management.MinimumBoundingGeometry(
        in_features=buffer,
        out_feature_class=mbp,
        geometry_type="ENVELOPE",
        group_option="NONE",
        group_field=None,
        mbg_fields_option="NO_MBG_FIELDS"
    )

    global mbp_project

    if veg_structure_crs != raster_crs:
        arcpy.Project_management(mbp, mbp_project, raster_crs)
    else:
        mbp_project = mbp


def create_NEON_plot_summary_stats():
    print("\nCalculate veg structure summary statistics (shrub and tree counts from CSV)...")
    arcpy.conversion.TableToTable(
        in_rows=apparent_individual_csv,
        out_path=intermediate_gdb,
        out_name=apparent_individual_table_name,
        where_clause="",
        config_keyword=""
    )

    arcpy.conversion.TableToTable(
        in_rows=apparent_individual_table,
        out_path=intermediate_gdb,
        out_name=shrubs_table_name,
        where_clause="growthForm LIKE '%shrub%'",
        config_keyword=""
    )

    arcpy.conversion.TableToTable(
        in_rows=apparent_individual_table,
        out_path=intermediate_gdb,
        out_name=trees_table_name,
        where_clause="growthForm LIKE '%tree%'",
        config_keyword=""
    )

    arcpy.analysis.Statistics(
        in_table=shrubs_table,
        out_table=summary_statistics_shrubs_i,
        statistics_fields="height MEAN",
        case_field="plotID;individualID",
        concatenation_separator=""
    )

    arcpy.analysis.Statistics(
        in_table=trees_table,
        out_table=summary_statistics_trees_i,
        statistics_fields="height MEAN",
        case_field="plotID;individualID",
        concatenation_separator=""
    )

    arcpy.AlterField_management(summary_statistics_shrubs_i, "MEAN_height", "shrub_height", "shrub_height")
    arcpy.AlterField_management(summary_statistics_trees_i, "MEAN_height", "tree_height", "tree_height")

    arcpy.analysis.Statistics(
        in_table=summary_statistics_shrubs_i,
        out_table=summary_statistics_shrubs,
        statistics_fields="shrub_height MEAN;shrub_height MEDIAN;shrub_height MIN;shrub_height MAX;shrub_height STD;shrub_height VARIANCE",
        case_field="plotID",
        concatenation_separator=""
    )

    arcpy.analysis.Statistics(
        in_table=summary_statistics_trees_i,
        out_table=summary_statistics_trees,
        statistics_fields="tree_height MEAN;tree_height MEDIAN;tree_height MIN;tree_height MAX;tree_height STD;tree_height VARIANCE",
        case_field="plotID",
        concatenation_separator=""
    )

    arcpy.AlterField_management(summary_statistics_shrubs, "FREQUENCY", "COUNT_shrubs", "COUNT_shrubs")
    arcpy.AlterField_management(summary_statistics_trees, "FREQUENCY", "COUNT_trees", "COUNT_trees")


def read_NEON_plot_summary_stats():
    print("\nRead veg structure summary statistics...")
    with arcpy.da.SearchCursor(summary_statistics_shrubs, "*") as sc:
        for row in sc:
            plot_id = row[1]
            if (plot_filters and plot_id in plot_filters) or (not plot_filters):
                count = row[2]
                mean = round(row[3], 2)
                median = round(row[4], 2)
                min = round(row[5], 2)
                max = round(row[6], 2)
                std = round(row[7], 2)
                var = round(row[8], 2)

                plot_dict[plot_id]["shrub_count"] = count
                plot_dict[plot_id]["shrub_h_mean"] = mean
                plot_dict[plot_id]["shrub_h_median"] = median
                plot_dict[plot_id]["shrub_h_min"] = min
                plot_dict[plot_id]["shrub_h_max"] = max
                plot_dict[plot_id]["shrub_h_std"] = std
                plot_dict[plot_id]["shrub_h_variance"] = var

    with arcpy.da.SearchCursor(summary_statistics_trees, "*") as sc:
        for row in sc:
            plot_id = row[1]
            if (plot_filters and plot_id in plot_filters) or (not plot_filters):
                count = row[2]
                mean = round(row[3], 2)
                median = round(row[4], 2)
                min = round(row[5], 2)
                max = round(row[6], 2)
                std = round(row[7], 2)
                var = round(row[8], 2)

                plot_dict[plot_id]["tree_count"] = count
                plot_dict[plot_id]["tree_h_mean"] = mean
                plot_dict[plot_id]["tree_h_median"] = median
                plot_dict[plot_id]["tree_h_min"] = min
                plot_dict[plot_id]["tree_h_max"] = max
                plot_dict[plot_id]["tree_h_std"] = std
                plot_dict[plot_id]["tree_h_variance"] = var


def extract_inversion_data():
    print("\nExtract cell values from inversion geotiffs...")

    rasters = arcpy.ListRasters()

    with arcpy.da.SearchCursor(mbp_project, ["SHAPE@", "plotID"]) as sc:
        for row in sc:
            extent = row[0].extent
            plot_id = row[1]
            if (plot_filters and plot_id in plot_filters) or (not plot_filters):
                #arcpy.MakeFeatureLayer_management(mbp_project, "mbp_project_select_" + plot_id, "plotID = '" + plot_id + "'")
                #arcpy.CopyFeatures_management(mbp_project_select, tmp_plot)
                #arcpy.env.extent=mbp_project
                tmp_plot = tmp_gdb + os.sep + "tmp_plot_" + plot_id
                expression = "plotID = '" + plot_id + "'"
                arcpy.analysis.Select(mbp_project, tmp_plot, expression)
                print("Plot ID: " + plot_id)
                temp_shrub_diff_dict = {}
                temp_tree_diff_dict = {}
                for raster in rasters:
                    if (raster_filters and any(x in raster for x in raster_filters)) or (not raster_filters):
                        #cells = arcpy.sa.ExtractByRectangle(raster, extent) # This gives cells with centroids outside the extent.
                        with arcpy.EnvManager(extent=tmp_plot):
                            cells = arcpy.sa.ExtractByMask(raster, tmp_plot)
                            save_raster = os.path.join(geotiffs_clip_dir, raster.split(".")[0] + "_clip_" + plot_id + ".tif")
                            cells.save(save_raster)
                        arr = arcpy.RasterToNumPyArray(save_raster, nodata_to_value=0)
                        sum = round(np.sum(arr), 2)
                        raster_basename = raster.split(".")[0]
                        plot_dict[plot_id][raster_basename] = sum

                        if "_Ns" in raster:
                            shrub_diff = round(plot_dict[plot_id]["shrub_count"] - sum, 2)
                            temp_shrub_diff_dict[raster_basename + "_diff"] = shrub_diff

                        if "_NT" in raster:
                            tree_diff = round(plot_dict[plot_id]["tree_count"] - sum, 2)
                            temp_tree_diff_dict[raster_basename + "_diff"] = tree_diff

                for k, v in temp_shrub_diff_dict.items():
                    plot_dict[plot_id][k] = v
                for k, v in temp_tree_diff_dict.items():
                    plot_dict[plot_id][k] = v


def write_to_csv():
    print("\nWrite to CSV...")
    df = pd.DataFrame.from_dict(plot_dict, orient="index")
    df.index.name = "plotID"
    df.to_csv(output_csv, index=True, header=True)

create_plot_fc()
create_NEON_plot_summary_stats()
read_NEON_plot_summary_stats()
extract_inversion_data()
write_to_csv()
