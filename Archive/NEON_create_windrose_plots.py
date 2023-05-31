from windrose import WindroseAxes
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import numpy as np
import csv
import pandas as pd
import os
import datetime

start_time = datetime.datetime.now()

#import arcpy
#arcpy.CheckOutExtension("GeoStats")
import datetime
#from arcpy.sa import *

wind_data_file = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Source\Wind\3d\NEON.D17.SOAP.DP1.00001.001.000.050.030.2DWSD_30min.2019-06.expanded.20221211T012100Z.csv"
output_folder = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\Charts\Windrose"


#df = pd.read_csv(wind_data_file, converters={'windSpeedMaximum': lambda x: 0 if x == "" else x, 'windDirMean': lambda x: 0 if x == "" else x}, usecols=['startDateTime','endDateTime','windSpeedMaximum','windDirMean'],keep_default_na=False)
df = pd.read_csv(wind_data_file, usecols=['startDateTime','endDateTime','windSpeedMaximum','windDirMean'], na_values=[""])

df = df.dropna()
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

plot_title = "SOAP"
plt.title(plot_title)

output_file = output_folder + os.sep  +"windrose.png"
plt.savefig(output_file)

end_time = datetime.datetime.now()
print "\nEnd Time: " + str(end_time)
print "Duration: " + str(end_time - start_time)