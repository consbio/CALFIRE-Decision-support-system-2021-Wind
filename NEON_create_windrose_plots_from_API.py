########################################################################################################################
# File name: NEON_create_windrose_plots_from_API.py
# Author: Mike Gough
# Date created: 05/31/2023
# Python Version: 2.7
# Description:
# Creates a windrose plot from 2D wind speed and direction data acquired through the NEON API:
# https://www.neonscience.org/resources/learning-hub/tutorials/neon-api-intro-requests-py
# More information on the data used in this script can be found at the following URL:
# https://data.neonscience.org/data-products/DP1.00001.001

#######################################################################################################################

from windrose import WindroseAxes
from matplotlib import pyplot as plt
import pandas as pd
import os
from datetime import datetime, date, time
import requests

# Input Parameters #####################################################################################################

NEON_site = "SOAP"

start_date = "2023-04-24"
end_date = "2023-04-27"

# Daytime (North winds, slow)
#start_hour = 03
#end_hour = 15

# Nighttime (South winds, fast)
start_hour = 3
end_hour = 15

interval = "30min"
sensor_position = "000.050"

output_dir = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\Charts\Windrose"

# NEON API Request #####################################################################################################

start_script = datetime.now()
print "\nStart Time: " + str(start_script)

SERVER = 'http://data.neonscience.org/api/v0/'
SITECODE = 'SOAP'

url = SERVER+'sites/'+SITECODE

PRODUCTCODE = 'DP1.00001.001'

start_date = datetime.strptime(start_date, "%Y-%m-%d")
end_date = datetime.strptime(end_date, "%Y-%m-%d")

year_str = str(start_date.year)
month_str = str(start_date.month).zfill(2)

data_request = requests.get(SERVER+'data/'+PRODUCTCODE+'/'+SITECODE+'/' + year_str + '-' + month_str)
data_json = data_request.json()

########################################################################################################################

for file in data_json['data']['files']: #Loop through keys of the data file dict
    if "csv" in file["name"].split(".")[-1]:
            if "basic" in file["name"]:
                if sensor_position == file["name"].split(".")[6] + "." + file["name"].split(".")[7]:
                    if interval == file["name"].split("_")[1].split(".")[0]:
                        print "\nNeon File Identified: "
                        for key, value in file.items():
                            print key, value
                        wind_data_file = file["url"]

df = pd.read_csv(wind_data_file, usecols=['startDateTime','endDateTime','windSpeedMaximum','windDirMean'], na_values=[""])

# Drop Null Values.
df = df.dropna()

# Filter the dataframe to only get dates between the user specified start and end date
df['startDate'] = pd.to_datetime(df['startDateTime'], format="%Y-%m-%d").dt.date
df['endDate'] = pd.to_datetime(df['endDateTime'], format="%Y-%m-%d").dt.date
df = df[(df['startDate'] >= start_date.date()) & (df['endDate'] <= end_date.date())]

# Filter the dataframe to only get hours between the user specified start and end hour
df['startTime'] = pd.to_datetime(df['startDateTime']).dt.time
df['endTime'] = pd.to_datetime(df['endDateTime']).dt.time
df = df[(df['startTime'] >= time(start_hour, 00)) & (df['endTime'] <= time(end_hour, 00))]

#start = df.startDateTime.tolist()
#end = df.endDateTime.tolist()
ws = map(float, df.windSpeedMaximum.tolist())
wd = map(float, df.windDirMean.tolist())

print("\nRows meeting requirements: \n" + df.to_string())

print "\nCreating windrose plot..."
ax = WindroseAxes.from_ax()
# It looks like the default angular labels are wrong!
# ax.set_xticklabels(['N', 'NW', 'W', 'SW', 'S', 'SE', 'E', 'NE'])
# Labels are off in windrose: https://github.com/python-windrose/windrose/issues/151
# ax.set_xticklabels([90, 45, 0, 315, 270, 225, 180, 135])
ax.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'])
ax.bar(wd, ws, normed=True, opening=0.8, edgecolor='white')
ax.set_legend(title='Wind Speed in m/s')

plot_title = "NEON FIELD SITE: " + NEON_site + "\n DATES: " + datetime.strftime(start_date, "%Y-%m-%d") + " - " + datetime.strftime(end_date, "%Y-%m-%d") + "  HOURS: " + str(start_hour).zfill(2) + "-" + str(end_hour).zfill(2) + " INTERVAL: " + interval + " SENSOR : " + sensor_position

ax.set_title(plot_title, pad=8)

output_file = output_dir + os.sep + "_".join(["NEON", NEON_site, "WINDROSE", datetime.strftime(start_date, "%Y%m%d"), datetime.strftime(end_date, "%Y%m%d"), str(start_hour).zfill(2), str(end_hour).zfill(2), interval, sensor_position]) + ".png"

plt.savefig(output_file)

# Open file with default OS application
os.system("start " + output_file)

end_script = datetime.now()
print "\nEnd Time: " + str(end_script)
print "Duration: " + str(end_script - start_script)