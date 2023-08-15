########################################################################################################################
# File name: NEON_create_windrose_plots_from_API.py
# Author: Mike Gough
# Date created: 08/14/2023
# Python Version: 3.x
# Description:
# Creates a csv and histogram from vegetation structure data (shrub height) acquired through the NEON API:
# https://www.neonscience.org/resources/learning-hub/tutorials/neon-api-intro-requests-py
# More information on the data used in this script can be found at the following URL:
# https://data.neonscience.org/data-products/DP1.10098.001
# Could potentially be modified in the future to get any statistics and plot data for any field.

#######################################################################################################################

import pandas as pd
import os
from datetime import datetime, date, time
import requests
import scipy
from scipy import stats
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

# Input Parameters #####################################################################################################

# NEON sites to process and date range.
NEON_sites = [
    ["SOAP", "2019-06-01", "2019-06-30"],
    ["TEAK", "2014-01-01", "2023-12-31"],
    ["SJER", "2014-01-01", "2023-12-31"],
]

# Vegetation structure.
PRODUCTCODE = 'DP1.10098.001'

# Output locations.
output_csv = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\CSV\stats\shrub_height_stats.csv"
output_histogram_dir = r"G:\CALFIRE_Decision_support_system_2021_mike_gough\Tasks\NEON\Data\Outputs\CSV\stats\histograms"

start_script = datetime.now()
print("\nStart Time: " + str(start_script))

stats_dict = {}

def get_shrub_stats():

    # NEON API Request #################################################################################################

    for NEON_site in NEON_sites:
        SITECODE = NEON_site[0]
        stats_dict[SITECODE] = {}
        print("\nSITE: " + SITECODE)
        start_date = datetime.strptime(NEON_site[1], "%Y-%m-%d")
        end_date = datetime.strptime(NEON_site[2], "%Y-%m-%d")

        start_date_str = start_date.strftime("%m/%d/%y")
        end_date_str = end_date.strftime("%m/%d/%y")

        print("Start Date: " + start_date_str)
        print("End Date: " + end_date_str)

        SERVER = 'http://data.neonscience.org/api/v0/'

        url = SERVER+'sites/'+SITECODE

        #Request the url
        site_request = requests.get(url)

        #Convert the request to Python JSON object
        site_json = site_request.json()
        #print(site_json)

        available_urls = []

        for product in site_json["data"]["dataProducts"]:
            # if a list item"s "dataProductCode" dict element equals the product code string
            if (product["dataProductCode"] == PRODUCTCODE):
                # print the available months
                print("Product: " + PRODUCTCODE)
                print("Available months: " + str(product["availableMonths"]))
                # print the available URLs

                available_months = []
                for month_and_url in zip(product['availableMonths'], product['availableDataUrls']):  # Loop through the list of months and URLs
                    available_month = datetime.strptime(month_and_url[0], "%Y-%m")
                    if available_month >= start_date and available_month <= end_date:
                        #print("URL: " + str(month_and_url[1]))
                        #print("Month: " + str(month_and_url[0]))
                        available_url = month_and_url[1]
                        available_urls.append(available_url)
                        available_months.append(month_and_url[0])

        print("Months between Start Date and End Date: " + str(available_months))

        df_list = []

        print("\nGetting CSV files for these months and reading into a Pandas Data Frame...")
        for request_url in available_urls:

            data_request = requests.get(request_url)
            data_json = data_request.json()
            #print(data_json)

            ####################################################################################################################

            for file_dict in data_json['data']['files']: #Loop through keys of the data file dict
                if "apparentindividual" in file_dict["name"]:
                    csv_url = file_dict["url"]
                    print(csv_url)
                    df_csv = pd.read_csv(csv_url, usecols=['uid', 'date', 'growthForm', 'height'], na_values=[""])
                    df_list.append(df_csv)

        print("\nCombining data frames...")
        df = pd.concat(df_list, ignore_index=True)

        # Drop Null Values.
        df = df.dropna()
        # print(df)

        #df = df[(df['growthForm'].str.contains("shrub"))]
        df['date'] = pd.to_datetime(df['date'], format="%Y-%m-%d").dt.date
        df = df[(df['growthForm'].str.contains("shrub")) & (df['date'] >= start_date.date()) & (df['date'] <= end_date.date())]

        #h = map(float, df.growthForm.tolist())
        height_list = df.height.tolist()

        if len(height_list) != 0:

            # Calculate statistics
            print("\nCalculating statistics...\n")

            #desc = scipy.stats.tvar(height_list, limits=None, inclusive=(True, True), axis=0, ddof=1)

            n = len(height_list)
            mn = min(height_list)
            mx = max(height_list)
            sm = sum(height_list)
            range = abs(mx-mn)
            sem = scipy.stats.sem(height_list, axis=0, ddof=1, nan_policy='propagate')
            med = scipy.stats.mstats.hdmedian(height_list, axis=-1, var=False)
            mode = scipy.stats.mode(height_list, axis=None).mode[0]
            mean = scipy.stats.tmean(height_list, limits=None, inclusive=(True, True), axis=None)
            kurt = scipy.stats.kurtosis(height_list, axis=0, fisher=True, bias=False, nan_policy='propagate')
            skew = scipy.stats.skew(height_list, axis=0, bias=False, nan_policy='propagate')
            var = scipy.stats.tvar(height_list, limits=None, inclusive=(True, True), axis=0, ddof=1)
            height_list_array = np.array(height_list)
            rms = np.sqrt(np.mean(height_list_array ** 2)) # https://stackoverflow.com/questions/40963659/root-mean-square-of-a-function-in-python
            std = scipy.stats.tstd(height_list, limits=None, inclusive=(True, True), axis=0, ddof=1)

            stats_dict[SITECODE]["start_date"] = start_date_str
            stats_dict[SITECODE]["end_date"] = end_date_str
            stats_dict[SITECODE]["months"] = str(available_months)
            stats_dict[SITECODE]["month_count"] = str(len(available_months))
            stats_dict[SITECODE]["mean"] = mean
            stats_dict[SITECODE]["sem"] = sem
            stats_dict[SITECODE]["med"] = med
            stats_dict[SITECODE]["mode"] = mode
            stats_dict[SITECODE]["std"] = std
            stats_dict[SITECODE]["var"] = var
            stats_dict[SITECODE]["kurt"] = kurt
            stats_dict[SITECODE]["skew"] = skew
            stats_dict[SITECODE]["range"] = range
            stats_dict[SITECODE]["min"] = mn
            stats_dict[SITECODE]["max"] = mx
            stats_dict[SITECODE]["sum"] = sm
            stats_dict[SITECODE]["count"] = n
            stats_dict[SITECODE]["rms"] = rms

            print("Mean: " + str(mean))
            print("Standard Error: " + str(sem))
            print("Median: " + str(med))
            print("Mode: " + str(mode))
            print("STD: " + str(std))
            print("Variance: " + str(var))
            print("Kurtosis: " + str(kurt))
            print("Skewness: " + str(skew))
            print("Range: " + str(range))
            print("Min : " + str(mn))
            print("Max : " + str(mx))
            print("Sum: " + str(sm))
            print("Count : " + str(n))
            print("RMS: " + str(rms))

            # Create Histogram
            print("\nCreating histogram...")

            df.hist(column='height', color='black')

            hist_title = SITECODE + " " + "Shrub Height"
            plt.suptitle(hist_title)
            plt.xlabel("Height (m)")
            plt.ylabel("Frequency")
            plt.legend(["Mean: " + str(round(mean,1))])

            hist_title = start_date_str + " - " + end_date_str + " (Months with Data: " + str(
                len(available_months)) + ")"
            plt.title(hist_title, fontsize=10)

            mean_label = str(round(mean,1))

            plt.axvline(mean, color='gray', linestyle='dashed', linewidth=1, label="Mean " + mean_label)
            #plt.legend(["Mean: " + str(round(mean,1))])
            plt.legend()

            output_histogram = os.path.join(output_histogram_dir, SITECODE + "_shrub_height_histogram.png")
            plt.savefig(output_histogram)
            # plt.show()

        else:
            print("No Shrubs")


get_shrub_stats()

stats_df = pd.DataFrame.from_dict(data=stats_dict, orient='index')
stats_df.index.name = 'Site'

print("Writing stats to CSV...")
stats_df.to_csv(output_csv, header=True)

print("Done.")

end_script = datetime.now()
print("\nEnd Time: " + str(end_script))
print("Duration: " + str(end_script - start_script))
