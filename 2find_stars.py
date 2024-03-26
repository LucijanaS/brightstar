"""
Created on: 28.02.2023
Created by: Lucijana Stanic

"""
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

from brightstar_input import lat_deg1, lon_deg1, height1, utc_offset, date_str
from brightstar_functions import dms_to_decimal

import csv
from datetime import datetime
from tqdm import tqdm
from itertools import zip_longest

import numpy as np

import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

# - * - coding: utf - 8 - * -

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

# Initialize a dictionary to store data from each column
data = {}

# Convert coordinates from deg to dec
lat_dec1 = dms_to_decimal(lat_deg1)
lon_dec1 = dms_to_decimal(lon_deg1)


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


# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

# Parse the input string into a datetime object
date_obj = datetime.strptime(date_str, '%Y-%m-%d')

# Convert the date object to a Julian date
date_JD = date_obj.toordinal() + 1721425 + .33333 - (
        1 / 24) * utc_offset  # added 1/3 since observations will most likely start at 8pm + offset of timezone

# Create a Time object from the observation time in Julian date
observation_time_utc = Time(date_JD, format='jd')

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------


# Extract stars which are visible on the night of observation
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
number_of_stars = 1000
print("Searching for stars visible in the night of", date_str)
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
        altaz_coords = equatorial_coord.transform_to(
            AltAz(obstime=time, location=EarthLocation(lat=lat_dec1, lon=lon_dec1, height=height1)))
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
    if low_altitude_count < len(altitudes) / 4:
        # Remove corresponding entries from the data dictionary for the star
        for key in data.keys():
            extracted_data[key].append(data[key][star_idx])

print("Out of the ", number_of_stars, " analysed, ", len(extracted_data[next(iter(data.keys()))]),
      " are visible throughout the night.")
list_1_print = input("Print list of those? (y or n) ")
if list_1_print == "y":
    print(extracted_data)

list_1_save = input("Save as .csv? (y or n) ")
if list_1_save == "y":
    output_file = 'stars_visible_' + date_str + '.csv'

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
