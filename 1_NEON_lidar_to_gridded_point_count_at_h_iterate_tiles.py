########################################################################################################################
# File name: 1_NEON_lidar_to_gridded_point_count_at_h_iterate.py
# Author: Mike Gough
# Date created: 06/01/2023
# Python Version: 3.x
# Description:
#
# For each tile a NEON study site, this script will generate a CSV containing the number of LiDAR point
# returns at different height above ground (HAG) intervals on a regular 2-dimensional grid (at a user specified
# resolution). For example, the field called "count_0m_1m" will store the number of point returns captured between 0-1m
# from the ground within each grid cell. Similarly, the field called "count_10m_11m will store the number of point
# returns captured between 10-11m from the ground within each grid cell.
#
# The HAG value for each LiDAR point return is determined by extracting the Digital Terrain Model (DTM)
# value to each point and then subtracting that value from the Z value of the point (HAG = Z - DTM).
#
# A list of the tiles to process comes from NEON_tiles_dataset (a reference shapefile provided by NEON).
#
# X,Y fields store the coordinate locations for the center of each grid cell.
# The output will snap to and match the coordinate system of the user defined snap_grid.
#
# Refer to the following URL for additional information on the NEON LiDAR data used:
# https://www.neonscience.org/data-collection/lidar
#
# Duration ~ 1.5 hours per NEON Tile.
########################################################################################################################

import arcpy
import os
from datetime import datetime
arcpy.env.overwriteOutput = True

# Input Parameters and Paths

version = "v8"  # Iterate all tiles

height_interval = 1
output_resolution = 10  # Units are m.
max_chm_offset = 30  # Max Canopy Height Offset allowed. LiDAR point returns greater than this distance above the CHM will be considered errors and will be deleted (units are m).

NEON_tiles_dataset = r"\\loxodonta\gis\Source_Data\environment\region\NEON_SITES\SOAP\Elevation_LiDAR\2021\07\NEON_lidar-elev\NEON.D17.SOAP.DP3.30024.001.2021-07.basic.20230601T181117Z.RELEASE-2023\NEON_D17_SOAP_DPQA_2021_merged_tiles.shp"
NEON_tiles_gdb = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Volume\NEON_Tiles\NEON_Tiles.gdb"
solar_insolation_index = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Source\Solar_Insolation\data\commondata\data0\insol_ind"
snap_grid = r"\\loxodonta\gis\Projects\CALFIRE_Decision_support_system_2021\Workspaces\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\CFO_Data_Processing\Data\Outputs\CFO_Green_Vegetation_By_County_Spring_2020\Fresno-County-California-Vegetation-GreenVegetation-2020-Spring-00010m.tif"
tmp_gdb = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Scratch\Scratch.gdb"
intermediate_gdb = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Volume\Volume_Iter_tiles.gdb"
intermediate_folder = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Volume"
output_fc_basename = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\Outputs.gdb\lidar_gridded_point_counts"
output_csv_folder = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\CSV"

start_script = datetime.now()
print("\nStart Time: " + str(start_script) + "\n")

print("### Getting a list of Tile IDs ###")

# Provide a list of NEON tiles to process, or create a list of all of them from the NEON_tiles_dataset
#tile_id_list = ["2021_SOAP_5_296000_4107000"]
tile_id_list = []

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
                NEON_DTM = r"\\loxodonta\gis\Source_Data\environment\region\NEON_SITES\SOAP\Elevation_LiDAR\2021\07\NEON_lidar-elev\NEON.D17.SOAP.DP3.30024.001.2021-07.basic.20230601T181117Z.RELEASE-2023\NEON_D17_SOAP_DP3_" + tile_id_selector + "_DTM.tif"
                NEON_CHM = r"\\loxodonta\gis\Source_Data\environment\region\NEON_SITES\SOAP\Ecosystem_structure\2019\06\NEON_struct-ecosystem\NEON.D17.SOAP.DP3.30015.001.2019-06.basic.20230524T172838Z.RELEASE-2023\NEON_D17_SOAP_DP3_" + tile_id_selector + "_CHM.tif"

                las_file_parse_date = datetime.strptime(NEON_lidar_laz_file.split("NEON.")[-1].split(".")[5] + '+0000', "%Y-%m%z")  # Specify Time Zone as GMT
                las_file_date = las_file_parse_date.strftime("%Y%m%d")
                extent_fc = os.path.join(NEON_tiles_gdb, "tile_" + tile_id)
                extent_name = extent_fc.split(os.sep)[-1]
                output_fc = "_".join([output_fc_basename, extent_name, las_file_date, version])
                output_csv_name = output_fc.split(os.sep)[-1] + ".csv"
                output_csv = os.path.join(output_csv_folder, output_csv_name)

                if not os.path.exists(NEON_lidar_laz_file):
                    print("No LiDAR file for this tile. Skipping.")
                elif not os.path.exists(NEON_DTM):
                    print("No DTM file")
                elif not os.path.exists(NEON_CHM):
                    print("No CHM file")
                elif os.path.exists(output_csv):
                    print("CSV file already exists. Skipping.")
                else:
                    tile_id_list.append(tile_id)

tile_count = len(tile_id_list)
print("\nTile Count: " + str(tile_count))

########################################################################################################################
count = 1

for tile_id in tile_id_list:

    print("\nTile: " + str(count) + "/" + str(tile_count) + "\n")
    print("\nTile ID: " + tile_id + "\n")

    tile_id_components = tile_id.split("_")
    tile_id_selector = "_".join([tile_id_components[-2], tile_id_components[-1]])

    # Paths to NEON Data (LiDAR, DTM, CHM) for this tile.
    NEON_lidar_laz_file = r"\\loxodonta\gis\Source_Data\environment\region\NEON_SITES\SOAP\Discrete_return_LiDAR_point_cloud\2019\06\NEON_lidar-point-cloud-line\NEON.D17.SOAP.DP1.30003.001.2019-06.basic.20230523T232633Z.RELEASE-2023\NEON_D17_SOAP_DP1_" + tile_id_selector + "_classified_point_cloud_colorized.laz"
    NEON_DTM = r"\\loxodonta\gis\Source_Data\environment\region\NEON_SITES\SOAP\Elevation_LiDAR\2021\07\NEON_lidar-elev\NEON.D17.SOAP.DP3.30024.001.2021-07.basic.20230601T181117Z.RELEASE-2023\NEON_D17_SOAP_DP3_" + tile_id_selector + "_DTM.tif"
    NEON_CHM = r"\\loxodonta\gis\Source_Data\environment\region\NEON_SITES\SOAP\Ecosystem_structure\2019\06\NEON_struct-ecosystem\NEON.D17.SOAP.DP3.30015.001.2019-06.basic.20230524T172838Z.RELEASE-2023\NEON_D17_SOAP_DP3_" + tile_id_selector + "_CHM.tif"

    if not os.path.exists(NEON_lidar_laz_file):
        print("No LiDAR file for this tile. Skipping.")
        continue

    print("### Creating Extent Feature Class for this Tile ###")

    # Test Extents
    # extent_fc = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Inputs\Extents\Extents.gdb\extent_1_tower_location"
    # extent_fc = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Inputs\Extents\Extents.gdb\SOAP_tile_298000_4100000"

    extent_fc = os.path.join(NEON_tiles_gdb, "tile_" + tile_id)
    expression = "TileID = '" + tile_id + "'"
    print(expression)
    arcpy.Select_analysis(NEON_tiles_dataset, extent_fc, expression)

    # Derived Parameters and Paths
    extent_name = extent_fc.split(os.sep)[-1]
    extent_fc_proj = os.path.join(intermediate_gdb, "extent_fc_proj_" + extent_name + "_" + version)

    lidar_folder = os.path.join(intermediate_folder, "Lidar")
    lidar_file_conversion_name = NEON_lidar_laz_file.split(os.sep)[-1].split(".")[0] + "_lidar_las_file_conversion.lasd"
    lidar_file_conversion = os.path.join(lidar_folder, lidar_file_conversion_name)

    lidar_clip_file = os.path.join(lidar_folder, lidar_file_conversion_name.split(".")[0] + "_clip" + ".lasd")
    input_las_file = os.path.join(lidar_folder, NEON_lidar_laz_file.split(os.sep)[-1].replace(".laz", ".las"))

    lidar_multipoint = os.path.join(intermediate_gdb, "lidar_multipoint_" + extent_name + "_" + version)

    lidar_points = os.path.join(intermediate_gdb, "lidar_points_" + extent_name + "_" + version)
    lidar_points_proj = os.path.join(intermediate_gdb, "lidar_points_proj_" + extent_name + "_" + version)

    extent_raster = os.path.join(intermediate_folder, "Extent_Rasters", "extent_raster_" + extent_name + "_" + version + ".tif")
    output_dtm = os.path.join(intermediate_gdb, "neon_dtm_clipped_" + extent_name + "_" + version)
    output_chm = os.path.join(intermediate_gdb, "neon_chm_clipped_" + extent_name + "_" + version)
    output_snap_grid = os.path.join(intermediate_gdb, "SNAP_GRID_clipped_" + extent_name + "_" + version)
    fishnet = os.path.join(intermediate_gdb, "fishnet_lidar_point_" + extent_name + "_" + version)

    tmp_points = os.path.join(tmp_gdb, "tmp_points")
    tmp_points_with_dtm = os.path.join(tmp_gdb, "tmp_points_with_dtm")
    tmp_points_with_dtm_and_chm = os.path.join(tmp_gdb, "tmp_points_with_dtm_and_chm")
    tmp_points_with_dtm_chm_and_insolation = os.path.join(tmp_gdb, "tmp_points_with_dtm_chm_and_insolation")

    las_file_parse_date = datetime.strptime(NEON_lidar_laz_file.split("NEON.")[-1].split(".")[5] + '+0000', "%Y-%m%z")  # Specify Time Zone as GMT
    las_file_date = las_file_parse_date.strftime("%Y%m%d")

    # Output Feature Class
    output_fc = "_".join([output_fc_basename, extent_name, las_file_date, version])
    output_csv_name = output_fc.split(os.sep)[-1] + ".csv"
    output_csv = os.path.join(output_csv_folder, output_csv_name)
    print(output_csv)
    if os.path.exists(output_csv):
        print("CSV file already exists. Skipping.")
        continue
    ########################################################################################################################

    #arcpy.env.extent = extent_fc
    arcpy.env.snapRaster = snap_grid

    extent_crs = arcpy.Describe(extent_fc).SpatialReference
    output_crs = arcpy.Describe(snap_grid).spatialReference


    def create_fishnet():
        """
        This function creates a fishnet from a raster version of the input extent snapped to the snap raster.
        """

        print("\n### Creating Fishnet ###\n")

        global extent_fc

        if extent_crs != output_crs:
            print("Projecting extent_fc to match CRS of snap_grid...")
            datum_extent_crs = extent_crs.gcs.name
            datum_output_crs = output_crs.gcs.name
            if datum_extent_crs != datum_output_crs:
                datum_transformation = "WGS_1984_(ITRF00)_To_NAD_1983"
            else:
                datum_transformation = ""

            arcpy.Project_management(
                in_dataset=extent_fc,
                out_dataset=extent_fc_proj,
                out_coor_system=output_crs,
                transform_method=datum_transformation,
                in_coor_system=extent_crs,
                preserve_shape="NO_PRESERVE_SHAPE", max_deviation="", vertical="NO_VERTICAL")

            extent_fc = extent_fc_proj

        #arcpy.env.snapRaster = snap_grid

        print("Converting extent_fc to raster and snapping to snap_grid...")
        with arcpy.EnvManager(snapRaster=snap_grid):
            arcpy.conversion.PolygonToRaster(extent_fc, "OBJECTID", extent_raster, "CELL_CENTER", "NONE", output_resolution, "BUILD")

        print("Creating Fishnet from extent_fc_raster...")
        desc = arcpy.Describe(extent_raster)
        arcpy.CreateFishnet_management(fishnet, str(desc.extent.lowerLeft),
                                       str(desc.extent.XMin) + " " + str(desc.extent.YMax), output_resolution, output_resolution, None, None,
                                       str(desc.extent.upperRight), "NO_LABELS", "#", "POLYGON")

        arcpy.management.DefineProjection(fishnet, output_crs)

        print("Removing Fishnet polys outside of extent_fc...")
        fishnet_layer = arcpy.MakeFeatureLayer_management(fishnet)
        arcpy.SelectLayerByLocation_management(fishnet_layer, "WITHIN", extent_fc, None, "NEW_SELECTION", "INVERT")
        arcpy.DeleteRows_management(fishnet_layer)

        print("Copying Fishnet to output feature class...")
        arcpy.CopyFeatures_management(fishnet, output_fc)

        print("Adding GRID ID to final output...")
        arcpy.AddField_management(output_fc, "grid_id", "LONG")
        arcpy.CalculateField_management(output_fc, "grid_id", "!OBJECTID!")


    def process_lidar_points():

        """
        This function converts the LiDAR Las file to a point feature class and adds the LiDAR Z value, DTM, and CHM value.
        Only points between the user defined starting height and the CHM + CHM offset variable will be included.
        Points below the DTM will be reassigned to height above ground = 0.
        """

        print("\n### Processing Lidar Points ###\n")
        print("LiDAR File Datetime: " + str(las_file_parse_date))

        global lidar_points

        print("Converting LAS file to an ArcGIS compatible LAS dataset file...")

        # Note that the LiDAR point processing all happens in the CRS native to the LiDAR and other NEON data.
        # This is desired for accuracy. Don't want to project the DTM, CHM, LiDAR. Extract in native CRS, then project.

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

        print("Converting LAS file...")

        arcpy.conversion.ConvertLas(
            NEON_lidar_laz_file,
            lidar_folder,
            "SAME_AS_INPUT", '', "NO_COMPRESSION", None,
            lidar_file_conversion,
            "NO_FILES", None)

        arcpy.ddd.ExtractLas(
            lidar_file_conversion,
            lidar_folder,
            "",
            fishnet, "PROCESS_EXTENT", '', "MAINTAIN_VLR", "REARRANGE_POINTS",
            "COMPUTE_STATS",
            lidar_clip_file,
            "SAME_AS_INPUT")

        print("Converting LAS file to Multipoints...")
        arcpy.ddd.LASToMultipoint(
            input_las_file,
            lidar_multipoint,
            1E-07, [], "ANY_RETURNS", None,
            #output_crs,
            "",
            "las", 1, "NO_RECURSION")

        print("Adding Z value to Multipoints...")
        arcpy.ddd.AddZInformation(lidar_multipoint, 'Z_max', 'NO_FILTER')

        print("Converting Multipoint to points...")  # Needed to extract by mask
        arcpy.management.FeatureToPoint(lidar_multipoint, lidar_points, "CENTROID")

        print("Extracting Raster values to points (DTM, CHM)...")
        # Note: these are extracted values for the LiDAR Points to calculate height above ground, max height, etc .
        # Since output resolution may be lower or different, doesn't make sense to use these in the final output.
        # Would need to take an average or do zonal statistics during post processing.
        arcpy.sa.ExtractMultiValuesToPoints(lidar_points, [[NEON_DTM, "dtm"], [NEON_CHM, "chm"]], "NONE")

        print("Calculating height of each point from ground (Z-Value - DTM)...")
        arcpy.AddField_management(lidar_points, "height_from_ground", "DOUBLE")
        with arcpy.da.UpdateCursor(lidar_points, ["Z_max", "dtm", "chm", "height_from_ground"]) as uc:
            for row in uc:
                if row[1] is not None and row[2] is not None: # Note that previously this was removing points where the CHM was 0 becuase used not row[1].
                    max_possible_height = row[2] + max_chm_offset  # Max height = CHM + max_chm_offset
                    height_from_ground = row[0] - row[1]  # Height from ground = Z value from LiDAR - DTM

                    # Delete LiDAR Point errors (points taller than the CHM)
                    # This is causing problems. Ground points are being deleted because they are above the CHM which is 0.
                    if height_from_ground > max_possible_height:
                        print("Deleting point with height = " + str(height_from_ground) + ". Max height from CHM is " + str(max_possible_height))
                        uc.deleteRow()
                    else:
                        # Assign all points below the ground a -1. This avoids having hundreds of negative fields.
                        if height_from_ground < 0:
                            height_from_ground = -1
                        row[3] = height_from_ground
                        uc.updateRow(row)
                else:
                    print("Deleting point without a DTM value")
                    uc.deleteRow()

        # If the lidar points are in a different CRS as the CFO grid (output crs), project them to match
        # for the spatial join in the next step.
        lidar_points_crs = arcpy.Describe(lidar_points).spatialReference
        if lidar_points_crs != output_crs:
            print("Projecting lidar points to match CRS of snap_grid...")
            datum_lidar_points = arcpy.Describe(lidar_points).spatialReference.gcs.name
            datum_output_crs = output_crs.gcs.name
            if datum_lidar_points != datum_output_crs:
                datum_transformation = "WGS_1984_(ITRF00)_To_NAD_1983"
            else:
                datum_transformation = ""

            arcpy.Project_management(
                in_dataset=lidar_points,
                out_dataset=lidar_points_proj,
                out_coor_system=output_crs,
                transform_method=datum_transformation,
                in_coor_system=lidar_points_crs,
                preserve_shape="NO_PRESERVE_SHAPE", max_deviation="", vertical="NO_VERTICAL")

            lidar_points = lidar_points_proj


    def count_point_returns():
        """
        This function counts the number of LiDAR point returns within regular vertical intervals measured from the ground
        (DTM). This information is stored in a Python dictionary, and added to a fishnet for the vertical slice using an
        update cursor.
        """

        print("\n### Counting Point Returns within Each Height Interval ###\n")

        all_h_values = [i[0] for i in arcpy.da.SearchCursor(lidar_points, "height_from_ground")]

        if len(all_h_values) == 0:
            print("No LiDAR points...")
            return 0

        h_min = int(min(all_h_values))
        h_max = int(round(max(all_h_values), 1))

        print("Min Height: " + str(h_min))
        print("Max Height: " + str(h_max) + "\n")

        h_start = h_min
        h_end = h_min + height_interval

        while h_end <= h_max + 1:

            print("Starting Height: " + str(h_start))
            print("Ending Height: " + str(h_end))

            print("Selecting liDAR points within this elevation range...")
            expression = "height_from_ground >= " + str(h_start) + " And height_from_ground < " + str(h_end)
            lidar_points_layer = arcpy.MakeFeatureLayer_management(lidar_points)
            arcpy.management.SelectLayerByAttribute(lidar_points_layer, "NEW_SELECTION", expression, None)

            print("Performing spatial join to get a count of the number LiDAR point returns within this height interval for each polygon...")
            output_spatial_join_fc = tmp_gdb + os.sep + "spatial_join_" + str(h_start).replace('-1', "neg") + "_" + str(h_end).replace("-1", "neg")
            arcpy.analysis.SpatialJoin(output_fc, lidar_points_layer, output_spatial_join_fc, "JOIN_ONE_TO_ONE", "KEEP_ALL", '' "CONTAINS", None, '')

            print("Creating a dictionary to store the point return count for each OBJECTID (OBJECTID[point_count])")
            count_dict = {}

            with arcpy.da.SearchCursor(output_spatial_join_fc, ["TARGET_FID", "Join_Count"]) as sc:
                for row in sc:
                    count_dict[row[0]] = row[1]

            #field_name = "count_" + str(h_start).replace("-1", "neg") + "m_" + str(h_end).replace("-1", "neg") + "m"
            if h_start == -1:
                field_name = "count_lt_0m"
            else:
                field_name = "count_" + str(h_start) + "m_" + str(h_end) + "m"
            print("Creating an update cursor to Add the point count to the following field in the output feature class: " + field_name + "\n")
            arcpy.AddField_management(output_fc, field_name, "LONG")
            with arcpy.da.UpdateCursor(output_fc, ["OBJECTID", field_name]) as uc:
                for row in uc:
                    row[1] = count_dict[row[0]]
                    uc.updateRow(row)

            h_start += 1
            h_end += 1


    def post_processing():
        """
        This function performs post-processing on the output feature class. Current operations include:
        1. Adding irregular height intervals (e.g., 1m-4m)
        2. Extracting raster values from a raster of interest (e.g., Solar Insolation).
        3. Adding x,y coordinates for the fishnet polygon centroids.
        """

        print("\n### Post Processing ###\n")

        all_h_values = [i[0] for i in arcpy.da.SearchCursor(lidar_points, "height_from_ground")]
        if len(all_h_values) == 0:
            print("No LiDAR points...")
            return 0

        h_max = int(round(max(all_h_values), 1))

        new_height_field = "count_1m_4m"
        print("Calculating new height field " + new_height_field)
        arcpy.AddField_management(output_fc, new_height_field, "LONG")
        height_fields = []
        for i in range(1, 4):
            height_fields.append("!count_" + str(i) + "m_" + str(i+1) + "m!")
        expression = " + ".join(height_fields)
        arcpy.CalculateField_management(output_fc, new_height_field, expression)

        new_height_field = "count_min_1m"
        print("Calculating new height field " + new_height_field)
        height_fields = []
        height_fields.append("!count_lt_0m!")
        height_fields.append("!count_0m_1m!")
        expression = " + ".join(height_fields)
        arcpy.AddField_management(output_fc, new_height_field, "LONG")
        arcpy.CalculateField_management(output_fc, new_height_field, expression)

        new_height_field = "count_5m_max"
        print("Calculating new height field " + new_height_field)
        height_fields = []
        for i in range(5, h_max):
            height_fields.append("!count_" + str(i) + "m_" + str(i+1) + "m!")
        expression = " + ".join(height_fields)
        arcpy.AddField_management(output_fc, new_height_field, "LONG")
        arcpy.CalculateField_management(output_fc, new_height_field, expression)


        print("Extracting Raster values to points (Insolation)...")
        arcpy.FeatureToPoint_management(output_fc, tmp_points)
        arcpy.sa.ExtractMultiValuesToPoints(tmp_points, [[solar_insolation_index, "insolation"]], "NONE")

        if False:
            #DTM is 1m, but final output may be Lower resolution (currently 10m). Probably need to do zonal stats here instead.
            print("Adding DTM to final output...")
            arcpy.sa.ExtractValuesToPoints(tmp_points, NEON_DTM, tmp_points_with_dtm)
            arcpy.AlterField_management(tmp_points_with_dtm, "RASTERVALU", "DTM_Extraction", "DTM_Extraction")

            print("Adding CHM to final output...")
            arcpy.sa.ExtractValuesToPoints(tmp_points_with_dtm, NEON_CHM, tmp_points_with_dtm_and_chm)
            arcpy.AlterField_management(tmp_points_with_dtm_and_chm, "RASTERVALU", "CHM_Extraction", "CHM_Extraction")


        arcpy.JoinField_management(output_fc, "grid_id", tmp_points, "grid_id", ["insolation"])

        print("Adding an X field...")
        arcpy.AddField_management(output_fc, "X", "Double")

        print("Adding a Y field...")
        arcpy.AddField_management(output_fc, "Y", "Double")

        print("Calculating X and Y Coordinates...")
        arcpy.management.CalculateGeometryAttributes(output_fc, "X CENTROID_X;Y CENTROID_Y", '', '', None, "SAME_AS_INPUT")


    def export_to_csv():
        arcpy.conversion.TableToTable(output_fc, output_csv_folder, output_csv_name)

    create_fishnet()
    process_lidar_points()
    count_point_returns()
    post_processing()
    export_to_csv()

    count += 1


end_script = datetime.now()
print("\nEnd Time: " + str(end_script))
print("Duration: " + str(end_script - start_script))
