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

# Input Parameters and Paths
version = "v1_10m_z11n"  # Match CRS to NEON SOAP (UTM Zone 10N)
version = "v2_10m_z10n"  # Match CRS to CFO data (UTM Zone 11N)
version = "v3_z10n_snap_trim"  # Set ENV snap raster to CFO, trim edges of fishnet where there is no data.
version = "v01_z10n_snap_trim"  # Set ENV snap raster to CFO, trim edges of fishnet where there is no data.
version = "v4_z10n"  # Full extent of the SOAP Tile containing the Tower.
version = "v5_z10n"  # Project input_points to match the CRS of the CFO data prior to spatial joins.

height_interval = 1
output_resolution = 10
max_chm_offset = 30  # Max Canopy Height Offset allowed. LiDAR point returns greater than this distance above the CHM will be considered errors and will be deleted (units are m).

#extent_fc = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Inputs\Extents\Extents.gdb\extent_1_tower_location"
#extent_fc = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Inputs\Extents\Extents.gdb\extent_1_tower_location_utm_zone_10n"
extent_fc = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Inputs\Extents\Extents.gdb\SOAP_tile_298000_4100000"
#extent_fc = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Inputs\Extents\Extents.gdb\SOAP_tile_298000_4100000_utm_zone_10n"

NEON_lidar_laz_file = r"\\loxodonta\gis\Source_Data\environment\region\NEON_SITES\SOAP\Discrete_return_LiDAR_point_cloud\2019\06\NEON_lidar-point-cloud-line\NEON.D17.SOAP.DP1.30003.001.2019-06.basic.20230523T232633Z.RELEASE-2023\NEON_D17_SOAP_DP1_298000_4100000_classified_point_cloud_colorized.laz"
NEON_DTM = r"\\loxodonta\gis\Source_Data\environment\region\NEON_SITES\SOAP\Elevation_LiDAR\2021\07\NEON_lidar-elev\NEON.D17.SOAP.DP3.30024.001.2021-07.basic.20230601T181117Z.RELEASE-2023\NEON_D17_SOAP_DP3_298000_4100000_DTM.tif"
NEON_CHM = r"\\loxodonta\gis\Source_Data\environment\region\NEON_SITES\SOAP\Ecosystem_structure\2019\06\NEON_struct-ecosystem\NEON.D17.SOAP.DP3.30015.001.2019-06.basic.20230524T172838Z.RELEASE-2023\NEON_D17_SOAP_DP3_298000_4100000_CHM.tif"
solar_insolation_index = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Source\Solar_Insolation\data\commondata\data0\insol_ind"

# Note: The output coordinate system will match the coordinate system of the snap_grid below
snap_grid = r"\\loxodonta\gis\Projects\CALFIRE_Decision_support_system_2021\Workspaces\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\CFO_Data_Processing\Data\Outputs\CFO_Green_Vegetation_By_County_Spring_2020\Fresno-County-California-Vegetation-GreenVegetation-2020-Spring-00010m.tif"

tmp_gdb = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Scratch\Scratch.gdb"
intermediate_gdb = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Volume\Volume.gdb"
intermediate_folder = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Volume"

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
fishnet_input_points_extent = os.path.join(intermediate_gdb, "fishnet_lidar_point_" + extent_name + "_" + version)

tmp_points = os.path.join(tmp_gdb, "tmp_points")
tmp_points_with_dtm = os.path.join(tmp_gdb, "tmp_points_with_dtm")
tmp_points_with_dtm_and_chm = os.path.join(tmp_gdb, "tmp_points_with_dtm_and_chm")
tmp_points_with_dtm_chm_and_insolation = os.path.join(tmp_gdb, "tmp_points_with_dtm_chm_and_insolation")

las_file_parse_date = datetime.strptime(NEON_lidar_laz_file.split("NEON.")[-1].split(".")[5] + '+0000', "%Y-%m%z")  # Specify Time Zone as GMT
print("\nLiDAR File Datetime: " + str(las_file_parse_date))
las_file_date = las_file_parse_date.strftime("%Y%m%d")

# Output Feature Class
output_fc = "_".join([r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\Outputs.gdb\lidar_gridded_point_counts_", extent_name, las_file_date, version])

########################################################################################################################

start_script = datetime.now()
print("\nStart Time: " + str(start_script))

arcpy.env.extent = extent_fc
arcpy.env.snapRaster = snap_grid

output_crs = arcpy.Describe(snap_grid).spatialReference
arcpy.env.outputCoordinateSystem = output_crs


def pre_processing():
    """
    This function creates a fishnet (from a raster version of the input extent snapped to the snap raster) and converts
    the LiDAR Las file to a point feature class containing the LiDAR Z value, DTM, and CHM extractions.
    Only points between the user defined starting height and the CHM + CHM offset variable will be included.
    Points below the DTM will be reassigned to height above ground = 0.
    """

    print("\nPreprocessing\n")

    global extent_fc, lidar_points

    extent_crs = arcpy.Describe(extent_fc).SpatialReference

    if extent_crs != output_crs:
        print("Projecting extent_fc to match CRS of snap_grid...")
        datum = extent_crs.gcs.name
        if datum == "GCS_North_American_1983":
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

    arcpy.env.snapRaster = snap_grid

    print("\nConverting Extent Feature Class to Extent Raster (Snapped to Snap Raster)...")
    with arcpy.EnvManager(snapRaster=snap_grid):
        arcpy.conversion.PolygonToRaster(extent_fc, "OBJECTID", extent_raster, "CELL_CENTER", "NONE", output_resolution, "BUILD")

    print("Creating Fishnet from Extent Raster...")
    desc = arcpy.Describe(extent_raster)
    arcpy.CreateFishnet_management(fishnet_input_points_extent, str(desc.extent.lowerLeft),
                                   str(desc.extent.XMin) + " " + str(desc.extent.YMax), output_resolution, output_resolution, None, None,
                                   str(desc.extent.upperRight), "NO_LABELS", "#", "POLYGON")

    print("Removing Fishnet polys outside of the Extent Raster...")
    fishnet_input_points_extent_layer = arcpy.MakeFeatureLayer_management(fishnet_input_points_extent)
    arcpy.SelectLayerByLocation_management(fishnet_input_points_extent_layer, "WITHIN", extent_fc, None, "NEW_SELECTION", "INVERT")
    arcpy.DeleteRows_management(fishnet_input_points_extent_layer)

    arcpy.management.DefineProjection(fishnet_input_points_extent, output_crs)

    print("Copying Fishnet to output Feature Class...")
    arcpy.CopyFeatures_management(fishnet_input_points_extent, output_fc)

    print("\nConverting LAS file to an ArcGIS compatible LAS dataset file...")

    #Note that the LAS file processing all happens in the CRS native to the NEON data.
    # This is desired for accuracy. Don't want to project the DTM, CHM, LiDAR. Extract in native CRS.

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

    print("Extracting LAS file to study area (and projecting LAS)...")
    arcpy.ddd.ExtractLas(
        lidar_file_conversion,
        lidar_folder,
        "",
        fishnet_input_points_extent_layer, "PROCESS_EXTENT", '', "MAINTAIN_VLR", "REARRANGE_POINTS",
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

    print("Calculating height of each point from ground (DTM)...")
    arcpy.AddField_management(lidar_points, "height_from_ground", "DOUBLE")
    with arcpy.da.UpdateCursor(lidar_points, ["Z_max", "dtm", "chm", "height_from_ground"]) as uc:
        for row in uc:
            if row[1]:
                max_possible_height = row[2] + max_chm_offset  # Max height = CHM + max_chm_offset
                height_from_ground = row[0] - row[1]  # Height from ground = Z value from LiDAR - DTM

                # Delete LiDAR Point errors (points taller than the CHM)
                # This is causing problems. Ground points are being deleted because they are above the CHM which is 0.
                if height_from_ground > max_possible_height:
                    print("Deleting point with height = " + str(height_from_ground) + ". Max height from CHM is " + str(max_possible_height))
                    uc.deleteRow()
                else:
                    if height_from_ground < 0:
                        height_from_ground = 0
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
        datum = extent_crs.gcs.name
        if datum == "GCS_North_American_1983":
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

    print("\nCounting Point Returns\n")

    all_h_values = [i[0] for i in arcpy.da.SearchCursor(lidar_points, "height_from_ground")]

    h_min = int(min(all_h_values))
    h_max = int(round(max(all_h_values), 1))

    print("Min Height: " + str(h_min))
    print("Max Height: " + str(h_max))

    h_start = h_min
    h_end = h_min + height_interval

    while h_end <= h_max + 1:

        print("Starting Height: " + str(h_start))
        print("Ending Height: " + str(h_end))

        print("Selecting liDAR points within this elevation range...")
        expression = "height_from_ground >= " + str(h_start) + " And height_from_ground < " + str(h_end)
        lidar_points_layer = arcpy.MakeFeatureLayer_management(lidar_points)
        arcpy.management.SelectLayerByAttribute(lidar_points_layer, "NEW_SELECTION", expression, None)

        print("Performing spatial join to get a count of the number LiDAR point returns within this height segment for each polygon...")
        output_spatial_join_fc = tmp_gdb + os.sep + "spatial_join_" + str(h_start).replace('-', "neg") + "_" + str(h_end).replace("-", "neg")
        arcpy.analysis.SpatialJoin(output_fc, lidar_points_layer, output_spatial_join_fc, "JOIN_ONE_TO_ONE", "KEEP_ALL", '' "CONTAINS", None, '')

        print("Creating a dictionary to store the return count for each OBJECTID (OBJECTID[point_count])")
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
    """
    This function performs post-processing on the output feature class. Current operations include:
    1. Adding irregular height intervals (e.g., 1m-4m)
    2. Extracting raster values from a raster of interest (e.g., Solar Insolation).
    3. Adding x,y coordinates for the fishnet polygon centroids.
    """

    print("\nPost Processing\n")

    all_h_values = [i[0] for i in arcpy.da.SearchCursor(lidar_points, "height_from_ground")]
    h_max = int(round(max(all_h_values), 1))

    new_height_field = "count_1m_4m"
    print("Calculating new height field " + new_height_field)
    arcpy.AddField_management(output_fc, new_height_field, "LONG")
    height_fields = []
    for i in range(1, 4):
        height_fields.append("!count_" + str(i) + "m_" + str(i+1) + "m!")
    expression = " + ".join(height_fields)
    arcpy.CalculateField_management(output_fc, new_height_field, expression)

    new_height_field = "count_5m_max"
    print("Calculating new height field " + new_height_field)
    height_fields = []
    for i in range(5, h_max):
        height_fields.append("!count_" + str(i) + "m_" + str(i+1) + "m!")
    expression = " + ".join(height_fields)
    arcpy.AddField_management(output_fc, new_height_field, "LONG")
    arcpy.CalculateField_management(output_fc, new_height_field, expression)

    print("Adding GRID ID to final output...")
    arcpy.AddField_management(output_fc, "GRID_ID", "LONG")
    arcpy.CalculateField_management(output_fc, "GRID_ID", "!OBJECTID!")

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


    arcpy.JoinField_management(output_fc, "GRID_ID", tmp_points, "GRID_ID", ["insolation"])

    print("Adding an X field...")
    arcpy.AddField_management(output_fc, "X", "Double")

    print("Adding a Y field...")
    arcpy.AddField_management(output_fc, "Y", "Double")

    print("Calculating X and Y Coordinates...")
    arcpy.management.CalculateGeometryAttributes(output_fc, "X CENTROID_X;Y CENTROID_Y", '', '', None, "SAME_AS_INPUT")

pre_processing()
count_point_returns()
post_processing()

end_script = datetime.now()
print("\nEnd Time: " + str(end_script))
print("Duration: " + str(end_script - start_script))
