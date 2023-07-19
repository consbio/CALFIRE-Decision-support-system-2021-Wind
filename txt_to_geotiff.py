########################################################################################################################
# File name: txt_to_geotiff.py
# Author: Mike Gough
# Date created: 07/19/2023
# Python Version: 3.x
# Description:
# Converts txt (csv) files into ascii rasters and then geotiffs.
# The txt (csv) files received don't have any spatial information in them. They are just txt (csv) files containing
# data values by row and column -- a subset of the CFO data extent (the row and column range was generated by the
# get_rows_and_cols.py script).
# This script opens up each text file, adds an ascii header, and replaces commas with spaces -- thereby making an ascii
# raster. This is then converted to a geotiff.
# All input files must be in the same SRS/CRS, and have the same extent and number of columns and rows.
# The number of rows and columns is determined automatically from the input file.
########################################################################################################################

import os
import shutil
import glob
from osgeo import gdal

input_dir = r"\\loxodonta\GIS\Source_Data\biota\region\NEON_Sites\SOAP\shrubs_and_trees_from_werne\little_tree_base_heights\b_lt_0"
ascii_dir = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Text_File_to_GeoTiff\ascii\b_lt_0"
geotiff_dir = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Text_File_to_GeoTiff\geotiff\b_lt_0"

input_dir = r"\\loxodonta\GIS\Source_Data\biota\region\NEON_Sites\SOAP\shrubs_and_trees_from_werne\little_tree_base_heights\b_lt_1"
ascii_dir = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Text_File_to_GeoTiff\ascii\b_lt_1"
geotiff_dir = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Intermediate\Text_File_to_GeoTiff\geotiff\b_lt_1"

# Spatial Reference System
srs = "EPSG:32610"  # WGS 84 / UTM zone 10N

# ASCII header params (ncols and nrows automatically determined from file)
# ncols = 1456
# nrows = 1425
xllcorner = 826260
yllcorner = 4099480
cellsize = 10
nodata_value = -32768

text_files = glob.glob(input_dir + "/*.txt")

column_count = None
row_count = None

for text_file in text_files:
    file_name = os.path.basename(text_file).split(".")[0]
    print("Input File: " + file_name)
    ascii_file = os.path.join(ascii_dir, file_name + ".asc")
    shutil.copyfile(text_file, ascii_file)
    with open(ascii_file, "r") as f:
        # Get the number of columns and rows from the first file.
        if not column_count and not row_count:
            lines = f.readlines()
            column_count = len(lines[0].split(","))
            row_count = len(lines)
            f.seek(0)  # start back at the beginning of the file.
            print("Column count: " + str(column_count))
            print("Row count: " + str(row_count))
        csv_text = f.read()
        ascii_text = csv_text.replace(",", " ")
    with open(ascii_file, "w") as f:
        f.write("ncols " + str(column_count) + "\n")
        f.write("nrows " + str(row_count) + "\n")
        f.write("xllcorner " + str(xllcorner) + "\n")
        f.write("yllcorner " + str(yllcorner) + "\n")
        f.write("cellsize " + str(cellsize) + "\n")
        f.write("nodata_value " + str(nodata_value) + "\n")
        f.write(ascii_text)
    geotiff_file = os.path.join(geotiff_dir, file_name + ".tif")
    ascii_ds = gdal.Open(ascii_file)
    gdal.Translate(geotiff_file, ascii_ds, outputSRS="EPSG:32610")
