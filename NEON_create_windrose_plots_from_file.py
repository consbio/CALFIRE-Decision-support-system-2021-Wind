from windrose import WindroseAxes
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import numpy as np
import csv
import pandas as pd
import os
from datetime import datetime, date, time
import requests

#import arcpy
#arcpy.CheckOutExtension("GeoStats")
#from arcpy.sa import *

start_date = datetime.strptime("2019-06-01", "%Y-%m-%d")
end_date = datetime.strptime("2019-06-20", "%Y-%m-%d")

start_hour = 12
end_hour = 03

start_script = datetime.now()

wind_data_file = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Source\Wind\3d\NEON.D17.SOAP.DP1.00001.001.000.050.030.2DWSD_30min.2019-06.expanded.20221211T012100Z.csv"
output_folder = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\Charts\Windrose"

#df = pd.read_csv(wind_data_file, converters={'windSpeedMaximum': lambda x: 0 if x == "" else x, 'windDirMean': lambda x: 0 if x == "" else x}, usecols=['startDateTime','endDateTime','windSpeedMaximum','windDirMean'],keep_default_na=False)
df = pd.read_csv(wind_data_file, usecols=['startDateTime','endDateTime','windSpeedMaximum','windDirMean'], na_values=[""])
df = df.dropna()

df['startDate'] = pd.to_datetime(df['startDateTime'], format="%Y-%m-%d").dt.date
df['endDate'] = pd.to_datetime(df['endDateTime'], format="%Y-%m-%d").dt.date

df['startTime'] = pd.to_datetime(df['startDateTime']).dt.time
df['endTime'] = pd.to_datetime(df['endDateTime']).dt.time

print(df['endDate'])
df = df[(df['startDate'] >= start_date.date()) & (df['endDate'] <= end_date.date())]
df = df[(df['startTime'] >= time(start_hour, 00)) & (df['endTime'] <= time(end_hour, 00))]

start = df.startDateTime.tolist()
end = df.endDateTime.tolist()
ws = map(float, df.windSpeedMaximum.tolist())
wd = map(float, df.windDirMean.tolist())

print start
print end
print ws
print wd

print "\nCreating windrose plot..."
# Create windrose plots.
# ws = np.random.random(500) * 6
# wd = np.random.random(500) * 360
#ws = np.array(wind_speed)
#wd = np.array(wind_dir_mean)
ax = WindroseAxes.from_ax()
# It looks like the default angular labels are wrong!
# ax.set_xticklabels(['N', 'NW', 'W', 'SW', 'S', 'SE', 'E', 'NE'])
# Labels are off in windrose: https://github.com/python-windrose/windrose/issues/151
# ax.set_xticklabels([90, 45, 0, 315, 270, 225, 180, 135])
ax.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'])
ax.bar(wd, ws, normed=True, opening=0.8, edgecolor='white')
ax.set_legend(title='Wind Speed in m/s')

#NEON_site = wind_data_file.split(os.sep)[-1].split(".")[2]
NEON_site = "SOAP"

plot_title = "NEON FIELD SITE: " + NEON_site + "\n DATES: " + datetime.strftime(start_date, "%Y-%m-%d") + "-" + datetime.strftime(end_date, "%Y-%m-%d") + "  HOURS: " + str(start_hour).zfill(2) + "-" + str(end_hour).zfill(2)
#plt.suptitle(plot_title)

#ttl = ax.title
#ttl.set_position([.5, 2])
ax.set_title(plot_title, pad=8)
ax.set_title(plot_title, pad=8)

output_file = output_folder + os.sep + "NEON_" + NEON_site + "_WINDROSE_" + datetime.strftime(start_date, "%Y%m%d") + "-" + datetime.strftime(end_date, "%Y%m%d") + "_" + str(start_hour).zfill(2) + "_" + str(end_hour).zfill(2) + ".png"

plt.savefig(output_file)

end_script = datetime.now()
print "\nEnd Time: " + str(end_script)
print "Duration: " + str(end_script - start_script)