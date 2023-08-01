########################################################################################################################
# File name: NEON_intensity_ratio.py
# Author: Mike Gough
# Date created: 07/12/2023
# Python Version: 3.x
# Description:
# Creates a ratio of ground return intensity to vegetation return intensity.
# Duration ~ 1.5 hours per NEON Tile.
########################################################################################################################

import arcpy
import os
from datetime import datetime
arcpy.env.overwriteOutput = True

# Input Parameters and Paths

version = "v1"  # Iterate all tiles

height_interval = 1
output_resolution = 10  # Units are m.

NEON_tiles_dataset = r"\\loxodonta\gis\Source_Data\environment\region\NEON_SITES\SOAP\Elevation_LiDAR\2021\07\NEON_lidar-elev\NEON.D17.SOAP.DP3.30024.001.2021-07.basic.20230601T181117Z.RELEASE-2023\NEON_D17_SOAP_DPQA_2021_merged_tiles.shp"
NEON_tiles_gdb = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Volume\NEON_Tiles\NEON_Tiles.gdb"


snap_grid = r"\\loxodonta\gis\Projects\CALFIRE_Decision_support_system_2021\Workspaces\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\CFO_Data_Processing\Data\Outputs\CFO_Green_Vegetation_By_County_Spring_2020\Fresno-County-California-Vegetation-GreenVegetation-2020-Spring-00010m.tif"
tmp_gdb = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Scratch\Scratch.gdb"
intermediate_gdb = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Intensity\Intensity.gdb"
intermediate_folder = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Intensity"

output_basename = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\Outputs.gdb\lidar_gridded_point_counts"

output_folder = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\Intensity"

output_crs = arcpy.Describe(snap_grid).spatialReference

start_script = datetime.now()
print("\nStart Time: " + str(start_script) + "\n")

print("### Getting a list of Tile IDs ###")

# Get a list of NEON tiles to process, or create a list of all of them from the NEON_tiles_dataset
tile_id_list = []
tile_id_list = ["2021_SOAP_5_300000_4096000"]

if not tile_id_list or len(tile_id_list) == 0:
    tile_id_list = []
    with arcpy.da.SearchCursor(NEON_tiles_dataset, ["TileID", "SHAPE@AREA"]) as sc:
        for row in sc:
            tile_id = row[0]
            tile_area = row[1]
            if tile_id != "merged_tiles" and tile_area > 999968:  # Should be 999968 to get all of the full tiles.
                # Paths to NEON Data (LiDAR, DTM, CHM) for this tile.
                tile_id_components = tile_id.split("_")
                tile_id_selector = "_".join([tile_id_components[-2], tile_id_components[-1]])
                NEON_lidar_laz_file = r"\\loxodonta\gis\Source_Data\environment\region\NEON_SITES\SOAP\Discrete_return_LiDAR_point_cloud\2019\06\NEON_lidar-point-cloud-line\NEON.D17.SOAP.DP1.30003.001.2019-06.basic.20230523T232633Z.RELEASE-2023\NEON_D17_SOAP_DP1_" + tile_id_selector + "_classified_point_cloud_colorized.laz"


                las_file_parse_date = datetime.strptime(NEON_lidar_laz_file.split("NEON.")[-1].split(".")[5] + '+0000', "%Y-%m%z")  # Specify Time Zone as GMT
                las_file_date = las_file_parse_date.strftime("%Y%m%d")

                base_name = "_".join([tile_id, las_file_date, version])
                intensity_ground_return = os.path.join(intermediate_folder, "intensity_ground_return_" + base_name + ".tif")
                intensity_veg_return = os.path.join(intermediate_folder, "intensity_veg_return_" + base_name + ".tif")
                intensity_ratio = os.path.join(output_folder, "intensity_ratio_ground_to_veg_" + base_name + ".tif")

                if not os.path.exists(NEON_lidar_laz_file):
                    print("No LiDAR file for this tile. Skipping.")
                elif os.path.exists(intensity_ratio):
                    print("Intensity Ratio already exists. Skipping.")
                else:
                    tile_id_list.append(tile_id)

tile_count = len(tile_id_list)
print("\nTile Count: " + str(tile_count))

########################################################################################################################

count = 1

for tile_id in tile_id_list:

    print("Create ratio of ground return to veg return intensity...")

    print("\nCount: " + str(count) + "/" + str(tile_count))
    print("\nTile ID: " + tile_id + "\n")

    tile_id_components = tile_id.split("_")
    tile_id_selector = "_".join([tile_id_components[-2], tile_id_components[-1]])

    # Paths to NEON LiDAR Data for this tile.
    NEON_lidar_laz_file = r"\\loxodonta\gis\Source_Data\environment\region\NEON_SITES\SOAP\Discrete_return_LiDAR_point_cloud\2019\06\NEON_lidar-point-cloud-line\NEON.D17.SOAP.DP1.30003.001.2019-06.basic.20230523T232633Z.RELEASE-2023\NEON_D17_SOAP_DP1_" + tile_id_selector + "_classified_point_cloud_colorized.laz"

    if not os.path.exists(NEON_lidar_laz_file):
        print("No LiDAR file for this tile. Skipping.")
        continue

    intensity_ground_returns = os.path.join(intermediate_folder, "intensity_ground_returns_" + tile_id + ".tif")
    intensity_veg_returns = os.path.join(intermediate_folder, "intensity_veg_returns_" + tile_id + ".tif")
    intensity_ratio_output = os.path.join(output_folder, "intensity_ratio_ground_to_veg_" + tile_id + ".tif")

    lidar_tmp_folder = os.path.join(intermediate_folder, "temp", "Lidar")

    lidar_file_lasd_conversion_name = NEON_lidar_laz_file.split(os.sep)[-1].split(".")[0] + "_lidar_las_file_conversion.lasd"
    lidar_file_lasd_conversion = os.path.join(lidar_tmp_folder, lidar_file_lasd_conversion_name)

    lidar_ground_returns_lasd_file_name = NEON_lidar_laz_file.split(os.sep)[-1].split(".")[0] + "_lidar_ground_returns_las_file_conversion.lasd"
    lidar_ground_returns_lasd_file = os.path.join(lidar_tmp_folder, lidar_ground_returns_lasd_file_name)

    lidar_veg_returns_lasd_file_name = NEON_lidar_laz_file.split(os.sep)[-1].split(".")[0] + "_lidar_veg_returns_las_file_conversion.lasd"
    lidar_veg_returns_lasd_file = os.path.join(lidar_tmp_folder, lidar_veg_returns_lasd_file_name)

    print("Convert LAS file to LASD file...")
    arcpy.conversion.ConvertLas(
        NEON_lidar_laz_file,
        lidar_tmp_folder,
        "SAME_AS_INPUT", '', "NO_COMPRESSION", None,
        lidar_file_lasd_conversion,
        "NO_FILES", None)

    #with arcpy.EnvManager(snapRaster=snap_grid, outputCoordinateSystem=output_crs, extent=lidar_file_lasd_conversion):
    with arcpy.EnvManager(snapRaster=snap_grid, outputCoordinateSystem=output_crs, extent=NEON_tiles_dataset):

        print("Make LASD layer consisting of only ground returns...")
        ground_returns = arcpy.management.MakeLasDatasetLayer(
            lidar_file_lasd_conversion, "ground_returns", "2", None,
            "INCLUDE_UNFLAGGED", "INCLUDE_SYNTHETIC", "INCLUDE_KEYPOINT", "EXCLUDE_WITHHELD", None, "INCLUDE_OVERLAP")

        input_crs = arcpy.Describe(lidar_file_lasd_conversion).spatialReference

        if input_crs != output_crs:
            datum_extent_crs = input_crs.gcs.name
            datum_output_crs = output_crs.gcs.name
            if datum_extent_crs != datum_output_crs:
                print("Setting geographic transformation...")
                arcpy.env.geographicTransformations = "WGS_1984_(ITRF00)_To_NAD_1983"


            print("Create raster of ground return intensity...")
            arcpy.conversion.LasDatasetToRaster(ground_returns, intensity_ground_returns, "INTENSITY",
                                                "BINNING AVERAGE LINEAR", "FLOAT", "CELLSIZE", output_resolution, 1)

            print("Make LASD layer consisting of only veg returns...")
            veg_returns = arcpy.management.MakeLasDatasetLayer(
                lidar_file_lasd_conversion, "veg_returns", "3;4;5", None,
                "INCLUDE_UNFLAGGED", "INCLUDE_SYNTHETIC", "INCLUDE_KEYPOINT", "EXCLUDE_WITHHELD", None, "INCLUDE_OVERLAP")

            print("Create raster of veg return intensity...")
            arcpy.conversion.LasDatasetToRaster(veg_returns, intensity_veg_returns, "INTENSITY",
                                                "BINNING AVERAGE LINEAR", "FLOAT", "CELLSIZE", output_resolution, 1)

            intensity_ratio = arcpy.sa.Divide(intensity_ground_returns, intensity_veg_returns)

            intensity_ratio.save(intensity_ratio_output)



