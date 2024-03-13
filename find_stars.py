"""
Created on: 28.02.2023
Created by: Lucijana Stanic

"""
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

import numpy as np
import astropy.units as u
import csv
from datetime import date, datetime
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import pytz
from timezonefinder import TimezoneFinder
from astropy.time import Time, TimeDelta
from tqdm import tqdm
from geopy.geocoders import Nominatim
from itertools import zip_longest



# - * - coding: utf - 8 - * -

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

# Initialize a dictionary to store data from each column
data = {}

# Open the CSV file
with open('1000stars_data.csv', newline='') as csvfile:
    reader = csv.reader(csvfile)
    header = next(reader)  # Read the header row

    # Create an empty list for each column
    for column_name in header:
        data[column_name] = []

    # Iterate over each row in the CSV file
    for row in reader:
        # Iterate over each item in the row and append it to the respective column list
        for idx, item in enumerate(row):
            column_name = header[idx]  # Get the corresponding column name
            data[column_name].append(item)



#date_str = input('Enter observation date in YYYY-MM-DD: ')
date_str = "2024-02-27"

#lat = input("Latitude (DD): ")
#lon = input("Longitude (DD): ")
#height = input("Height (m above sea level, leave empty for no height): ")
lat=47.3769*u.deg
lon=8.5417*u.deg
height=0*u.m

if height is None:
    height = 0*u.m

# Parse the input string into a datetime object
date_obj = datetime.strptime(date_str, '%Y-%m-%d')

# Convert the date object to a Julian date
utc_offset = 1
date_JD = date_obj.toordinal() + 1721425 + .33333 - (1/24)*utc_offset # added 1/3 such since observations will most likely start at 8pm + offset of timezone

# Create a Time object from the observation time in Julian date
observation_time_utc = Time(date_JD, format='jd')
print('Enter observation site coordinates ')

RA = []
for star in data["RA_decimal"]:
    RA.append(star)

Dec = []
for star in data["Dec_decimal"]:
    Dec.append(star)

equatorial_coords = []
for i in range(len(data["RA_decimal"])):
    coord = SkyCoord(RA[i], Dec[i], unit=(u.hourangle, u.deg), frame='icrs')
    equatorial_coords.append(coord)

# Initialize lists to store altitude and azimuth values for each star
altitudes_per_star = []
azimuths_per_star = []

# Iterate over each equatorial coordinate
number_of_stars = 50
for equatorial_coord in tqdm(equatorial_coords[:number_of_stars]):
    # Initialize lists to store altitude and azimuth values for the current star
    altitudes = []
    azimuths = []

    # Define the time span
    hours_before = 0
    hours_after = 12
    start_time = observation_time_utc - TimeDelta(hours_before * u.hour)
    end_time = observation_time_utc + TimeDelta(hours_after * u.hour)

    # Calculate altitude and azimuth for each time point
    times = start_time + (end_time - start_time) * np.linspace(0, 1, 97)[:, None]
    for time in times:
        altaz_coords = equatorial_coord.transform_to(AltAz(obstime=time, location=EarthLocation(lat=lat, lon=lon, height=height)))
        altitude = altaz_coords.alt
        azimuth = altaz_coords.az
        altitudes.append(altitude)
        azimuths.append(azimuth)

    # Append altitude and azimuth lists for the current star to the main lists
    altitudes_per_star.append(altitudes)
    azimuths_per_star.append(azimuths)

# Define a threshold for altitude entries
altitude_threshold = 10  # Minimum altitude value during one night

# Iterate over each star's altitude data
extracted_data = {}
for key in data.keys():
    extracted_data[key] = []
for star_idx, altitudes in tqdm(enumerate(altitudes_per_star)):
    # Count the number of altitude entries less than the threshold
    low_altitude_count = sum(altitude.value < altitude_threshold for altitude in altitudes)

    # Check if more than 1/4 of the entries have altitudes less than the threshold
    if low_altitude_count > len(altitudes) / 4:
        # Remove corresponding entries for the star from altitudes_per_star and azimuths_per_star
        del altitudes_per_star[star_idx]
        del azimuths_per_star[star_idx]

        # Remove corresponding entries from the data dictionary for the star
        for key in data.keys():
            extracted_data[key].append(data[key][star_idx])

print("Out of the ", number_of_stars," analysed, ", len(extracted_data[next(iter(data.keys()))])," are visible throughout the night.")
list_1_print = input("Print list of those? (y or n) ")
if list_1_print == "y":
    print(extracted_data)

list_1_save = input("Save as .csv? (y or n) ")
if list_1_save == "y":
    output_file = 'stars_visible_'+date_str+'.csv'

    # Write the extracted data to the CSV file
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Write the header row using the keys of the extracted_data dictionary
        writer.writerow(extracted_data.keys())

        # Write the data rows
        # Determine the number of rows based on the length of one of the lists in extracted_data
        num_rows = len(next(iter(extracted_data.values())))

        for i in range(num_rows):
            # Get the data for each column and write it to the CSV file
            row_data = [extracted_data[key][i] for key in extracted_data.keys()]
            writer.writerow(row_data)
    print("Saved as", output_file)



magnitude = "n"
"""
magnitude = input("Do you want to specify the apparent magnitude (V) of the star? ")
if magnitude == "y":
    magnitude_threshold = float(input("Apparent magnitude threshold <= "))

    # Initialize extracted_data dictionary
    extracted_data2 = {}
    for key in data.keys():
        extracted_data2[key] = []

    # Loop over each star's Vmag and check against the threshold
    for star, vmag in enumerate(extracted_data["Vmag"]):
        if float(vmag) <= magnitude_threshold:
            # Append corresponding entries for the star to extracted_data
            for key in extracted_data.keys():
                extracted_data2[key].append(data[key][star])

    print("Out of the ", number_of_stars," analysed, ", len(extracted_data2[next(iter(data.keys()))])," are visible throughout the night and have and apparent magnitude <= ", magnitude_threshold)
    list_2_print = input("Print list of those? (y or n) ")
    if list_2_print == "y":
        print(extracted_data2)

    list_2_save = input("Save as .csv? (y or n) ")
    if list_2_save == "y":
        output_file = 'stars_visible_bright_'+date_str+'.csv'

        # Write the extracted data to the CSV file
        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)

            # Write the header row using the keys of the extracted_data dictionary
            writer.writerow(extracted_data2.keys())

            # Write the data rows
            # Determine the number of rows based on the length of one of the lists in extracted_data
            num_rows = len(next(iter(extracted_data2.values())))

            for i in range(num_rows):
                # Get the data for each column and write it to the CSV file
                row_data = [extracted_data2[key][i] for key in extracted_data2.keys()]
                writer.writerow(row_data)
        print("Saved as", output_file)
    extracted_data = extracted_data2
"""



diameter_choice = input("Do you want to specify the diameter in mas of the star? ")

if diameter_choice == "y":
    diameter_threshold = float(input("Diameter in mas >= "))

    # Initialize extracted_data dictionary
    extracted_data3 = {}
    for key in data.keys():
        extracted_data3[key] = []

    # Loop over each star's Vmag and check against the threshold
    for idx, (diam_V, diam_B, diam_U) in enumerate(zip_longest(extracted_data["Diameter_V"], extracted_data["Diameter_B"],
                                                    extracted_data["Diameter_U"])):
        # Convert diameters to float, handling empty strings
        diam_V = float(diam_V) if diam_V else None
        diam_B = float(diam_B) if diam_B else None
        diam_U = float(diam_U) if diam_U else None

        range_num = 0.1

        # Check if any of the diameters are within the threshold range_num
        if (diam_V is not None and abs(diam_V - diameter_threshold) <= range_num) or \
                (diam_B is not None and abs(diam_B - diameter_threshold) <= range_num) or \
                (diam_U is not None and abs(diam_U - diameter_threshold) <= range_num):
            # Append corresponding entries for the star to extracted_data
            for key in data.keys():
                extracted_data3[key].append(extracted_data[key][idx])


    if magnitude == "y":
        print("Out of the ", number_of_stars," analysed, ", len(extracted_data3[next(iter(data.keys()))])," are visible throughout the night and have and apparent magnitude <= ", magnitude_threshold, "and a diameter of around ", diameter_threshold, " mas.")
    else:
        print("Out of the ", number_of_stars, " analysed, ", len(extracted_data3[next(iter(data.keys()))])," are visible throughout the night and have a diameter of around ", diameter_threshold, " mas.")
    list_3_print = input("Print list of those? (y or n) ")
    if list_3_print == "y":
        print(extracted_data3)

    list_3_save = input("Save as .csv? (y or n) ")
    if list_3_save == "y":
        if magnitude == "y":
            output_file = 'stars_visible_bright_big_'+date_str+'.csv'
        else:
            output_file = 'stars_visible_big_'+date_str+'.csv'

        # Write the extracted data to the CSV file
        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)

            # Write the header row using the keys of the extracted_data dictionary
            writer.writerow(extracted_data3.keys())

            # Write the data rows
            # Determine the number of rows based on the length of one of the lists in extracted_data
            num_rows = len(next(iter(extracted_data3.values())))

            for i in range(num_rows):
                # Get the data for each column and write it to the CSV file
                row_data = [extracted_data3[key][i] for key in extracted_data3.keys()]
                writer.writerow(row_data)
        print("Saved as", output_file)
