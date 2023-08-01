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
raster_resolution = 10
raster_columns = 26577
raster_rows = 35595

# CFO (SSN Study Area) raster extent
raster_left = 688080
raster_right = 953850
raster_top = 4275130
raster_bottom = 3919180

# NEON SOAP extent in UTM Zone 10N
boundary_left = 826265.190808
boundary_right = 840810.153080
boundary_top = 4113725.003810
boundary_bottom = 4099487.636019

column_start = int(math.floor((boundary_left - raster_left)/raster_resolution))
column_end = int(math.ceil((boundary_right - raster_left)/raster_resolution))

row_start = int(math.floor((raster_top - boundary_top)/raster_resolution))
row_end = int(math.ceil((raster_top - boundary_bottom)/raster_resolution))

print("Column Range: " + str(column_start) + " - " + str(column_end))
print("Row Range: " + str(row_start) + " - " + str(row_end))

#print("Column Start: " + str(column_start))

xllcorner = raster_left + (column_start * 10)
yllcorner = raster_top - (row_end * 10)

print("xllcorner: " + str(xllcorner))
print("yllcorner: " + str(yllcorner))

