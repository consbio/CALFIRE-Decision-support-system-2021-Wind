########################################################################################################################
# File name: NEON_CHM_stats_in_CFO_grid.py
# Author: Mike Gough
# Date created: 11/03/2023
# Python Version: 3.x
# Description:
# This script mosaics the NEON CHM data (Canopy Height) and summarized the data to the CFO grid using zonal statistics.
# Zonal statistics are performed in the CRS of the NEON CHM data (UTM 11N), and then joined to the CFO fishnet
# which is in UTM 10N. This was done in order to avoid projecting and resampling the CHM data.
########################################################################################################################

import arcpy
import glob
import os

arcpy.env.overwriteOutput = True

# Original NEON CHM data Tiles (Count = 225)
source_dir = r"\\loxodonta\gis\Source_Data\environment\region\NEON_SITES\SOAP\Ecosystem_structure\2019\06\NEON_struct-ecosystem\NEON.D17.SOAP.DP3.30015.001.2019-06.basic.20230524T172838Z.RELEASE-2023"

# Mosaic raster output:
mosaic_dir = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Inputs\Canopy_Height_Model\2019\06\tif"
CHM_mosaic_name = "NEON_D17_SOAP_DP3_CHM_MOSAIC_2019-06.tif"

intermediate_dir = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Canopy_Height_Model"
intermediate_gdb = os.path.join(intermediate_dir, "Canopy_Height_Model.gdb")

output_dir = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\Canopy_Height_Model"

output_snap_raster = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\CFO_Calcs_From_LiDAR\GeoTiff\v1\neon_soap_lidar_to_canopy_cover_v1.tif"
UTM_11N_proj_file = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\WGS 1984 UTM Zone 11N.prj"

CFO_zone_fishnet = os.path.join(intermediate_dir, intermediate_gdb, "CFO_NEON_SOAP_Fishnet_Copy")
CFO_zone_fishnet_11N = os.path.join(intermediate_dir, intermediate_gdb, "CFO_NEON_SOAP_Fishnet_11N")

CHM_mosaic_projected_name = "NEON_D17_SOAP_DP3_CHM_MOSAIC_2019-06_proj.tif"
CHM_mosaic_projected_raster = os.path.join(intermediate_dir, CHM_mosaic_projected_name)

CHM_output_raster_max = os.path.join(output_dir, "neon_soap_chm_mosaic_2019-06_cfo_zonal_max.tif")
CHM_output_raster_std = os.path.join(output_dir, "neon_soap_chm_mosaic_2019-06_cfo_zonal_std.tif")

arcpy.env.workspace = mosaic_dir
CHM_rasters = glob.glob(source_dir + os.sep + "*CHM.tif")

CFO_zone_fishnet_11N_raster = os.path.join(intermediate_dir, "CFO_NEON_SOAP_Fishnet_11N_Raster_1m_Snapped_CHM.tif")

CHM_zonal_stats_table = os.path.join(intermediate_gdb, "CHM_Zonal_Stats_Table")


def mosaic():
    print("Mosaicing rasters...")
    print("Count: " + str(len(CHM_rasters))) # 225 rasters for NEON SOAP
    arcpy.MosaicToNewRaster_management(CHM_rasters, mosaic_dir, CHM_mosaic_name, "", "32_BIT_FLOAT", 1, 1)


def create_fishnet_zones():
    # Zone fishnet extent parameters set to match canopy_cover raster created from LiDAR data.
    # G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\CFO_Calcs_From_LiDAR\GeoTiff\v1\neon_soap_lidar_to_canopy_cover_v1.tif
    print("Creating fishnet using canopy_cover created from lidar data that's on the CFO grid...")
    arcpy.management.CreateFishnet(
        out_feature_class=CFO_zone_fishnet,
        origin_coord="826650.00012207 4099790.00012207",
        y_axis_coord="826650.00012207 4099800.00012207",
        cell_width=10,
        cell_height=10,
        number_rows=None,
        number_columns=None,
        corner_coord="839810.00012207 4113580.00012207",
        labels="NO_LABELS",
        template='826650.00012207 4099790.00012207 839810.00012207 4113580.00012207 PROJCS["WGS_1984_UTM_Zone_10N",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",-123.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]',
        geometry_type="POLYGON"
    )

    print("Adding ORIG_OID field...")
    arcpy.AddField_management(CFO_zone_fishnet, "ORIG_OID", "LONG")
    arcpy.CalculateField_management(CFO_zone_fishnet, "ORIG_OID", "!OID!")


def project_fishnet():
    print("Projecting fishnet to UTM 11N to match NEON CHM data...")
    arcpy.Project_management(CFO_zone_fishnet, CFO_zone_fishnet_11N, UTM_11N_proj_file)


def convert_fishnet_to_raster():
    print("Converting projected fishnet to raster. Snapping to NEON CHM data at same resolution 1m.")
    with arcpy.EnvManager(snapRaster=CHM_mosaic_name):
        arcpy.conversion.FeatureToRaster(
            in_features=CFO_zone_fishnet_11N,
            field="ORIG_OID",
            out_raster=CFO_zone_fishnet_11N_raster,
            cell_size=1
        )


def calc_zonal_stats_as_table():
    print("Calculating zonal statistics in 11N (original CHM CRS)...")
    arcpy.sa.ZonalStatisticsAsTable(CFO_zone_fishnet_11N_raster, "Value", CHM_mosaic_name, CHM_zonal_stats_table,  "DATA", "ALL")


def join_zonal_stats():
    print("Joining zonal statistics to original CFO fishnet in 10N...")
    arcpy.JoinField_management(CFO_zone_fishnet, "ORIG_OID", CHM_zonal_stats_table, "Value")


def zonal_stats_to_raster():
    print("Converting zonal statistics in original CFO fishnet to raster...")
    with arcpy.EnvManager(snapRaster=output_snap_raster, extent=output_snap_raster):
        arcpy.conversion.FeatureToRaster(
            in_features=CFO_zone_fishnet,
            field="MAX",
            out_raster=CHM_output_raster_max,
            cell_size=1
        )
        arcpy.conversion.FeatureToRaster(
            in_features=CFO_zone_fishnet,
            field="STD",
            out_raster=CHM_output_raster_std,
            cell_size=1
        )


mosaic()
create_fishnet_zones()
project_fishnet()
convert_fishnet_to_raster()
calc_zonal_stats_as_table()
join_zonal_stats()
zonal_stats_to_raster()




