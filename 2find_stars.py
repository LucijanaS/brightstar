"""
Created on: 28.02.2023
Created by: Lucijana Stanic

"""
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

from brightstar_input import lat_deg1, lon_deg1, height1, utc_offset, date_str, n_brightest_stars
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
    reader = csv.DictReader(csvfile)
    header = reader.fieldnames

    # Create an empty list for each column
    for column_name in header:
        data[column_name] = []

    # Iterate over each row in the CSV file
    for row in reader:
        for column_name in header:
            data[column_name].append(row[column_name])

# Parse the input string into a datetime object
date_obj = datetime.strptime(date_str, '%Y-%m-%d')
date_JD = date_obj.toordinal() + 1721425 + .33333 - (1 / 24) * utc_offset
observation_time_utc = Time(date_JD, format='jd')

# Extract stars which are visible on the night of observation
equatorial_coords = [
    SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='icrs')
    for ra, dec in zip(map(float, data["RA_decimal"]), map(float, data["Dec_decimal"]))
]

altitudes_per_star = []
azimuths_per_star = []

# Define the time span
hours_before = 0
hours_after = 12
start_time = observation_time_utc - TimeDelta(hours_before * u.hour)
end_time = observation_time_utc + TimeDelta(hours_after * u.hour)

for equatorial_coord in tqdm(equatorial_coords[:10]):
    altitudes, azimuths = [], []
    times = start_time + (end_time - start_time) * np.linspace(0, 1, 97)[:, None]
    for time in times:
        altaz_coords = equatorial_coord.transform_to(AltAz(obstime=time, location=EarthLocation(lat=lat_dec1, lon=lon_dec1, height=height1)))
        altitudes.append(altaz_coords.alt)
        azimuths.append(altaz_coords.az)
    altitudes_per_star.append(altitudes)
    azimuths_per_star.append(azimuths)

# Filter stars based on altitude threshold
altitude_threshold = 10
extracted_data = {key: [] for key in data.keys()}
for star_idx, altitudes in tqdm(enumerate(altitudes_per_star)):
    low_altitude_count = sum(altitude.value < altitude_threshold for altitude in altitudes)
    if low_altitude_count < len(altitudes) / 4:
        for key in data.keys():
            extracted_data[key].append(data[key][star_idx])

# Save extracted data to CSV
print("Out of the ", n_brightest_stars, " analysed, ", len(extracted_data[next(iter(data.keys()))]),
      " are visible throughout the night.")
list_1_print = input("Print list of those? (y or n) ")
if list_1_print.lower() == "y":
    print(extracted_data)

list_1_save = input("Save as .csv? (y or n) ")
if list_1_save.lower() == "y":
    output_file = 'stars_visible_' + date_str + '.csv'
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=extracted_data.keys())
        writer.writeheader()
        writer.writerows({key: value[i] for key, value in extracted_data.items()} for i in range(len(altitudes_per_star)))

# Ask for user input
magnitude = input("Do you want to specify the apparent magnitude (V) of the star? ").lower()
if magnitude == "y":
    magnitude_threshold = float(input("Apparent magnitude threshold <= "))

    extracted_data2 = {key: [] for key in data.keys()}
    for star, vmag in zip_longest(extracted_data["Vmag"], extracted_data["Vmag"]):
        if vmag is not None and float(vmag) <= magnitude_threshold:
            for key in extracted_data.keys():
                extracted_data2[key].append(data[key][star_idx])

    print("Out of the ", n_brightest_stars, " analysed, ", len(extracted_data2[next(iter(data.keys()))]), "are visible "
                                                                                                        "throughout the "
                                                                                                        "night and have "
                                                                                                        "and apparent "
                                                                                                        "magnitude <= "
                                                                                                        "",
          magnitude_threshold)
    list_2_print = input("Print list of those? (y or n) ").lower()
    if list_2_print == "y":
        print(extracted_data2)

    list_2_save = input("Save as .csv? (y or n) ").lower()
    if list_2_save == "y":
        output_file = 'stars_visible_bright_' + date_str + '.csv'
        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=extracted_data2.keys())
            writer.writeheader()
            writer.writerows(extracted_data2)

    extracted_data = extracted_data2

# Ask for user input for diameter choice
diameter_choice = input("Do you want to specify the diameter in mas of the star? ").lower()

if diameter_choice == "y":
    diameter_threshold = float(input("Diameter in mas >= "))

    extracted_data3 = {key: [] for key in data.keys()}
    for idx, (diam_V, diam_B, diam_U) in enumerate(zip_longest(extracted_data["Diameter_V"], extracted_data["Diameter_B"], extracted_data["Diameter_U"])):
        diam_V = float(diam_V) if diam_V else None
        diam_B = float(diam_B) if diam_B else None
        diam_U = float(diam_U) if diam_U else None
        range_num = 0.1

        if (diam_V is not None and abs(diam_V - diameter_threshold) <= range_num) or \
                (diam_B is not None and abs(diam_B - diameter_threshold) <= range_num) or \
                (diam_U is not None and abs(diam_U - diameter_threshold) <= range_num):
            for key in extracted_data.keys():
                extracted_data3[key].append(extracted_data[key][idx])

    if magnitude == "y":
        print("Out of the ", n_brightest_stars, " analysed, ", len(extracted_data3[next(iter(data.keys()))]),
              " are visible throughout the night and have and apparent magnitude <= ", magnitude_threshold,
              "and a diameter of around ", diameter_threshold, " mas.")
    else:
        print("Out of the ", n_brightest_stars, " analysed, ", len(extracted_data3[next(iter(data.keys()))]),
              " are visible throughout the night and have a diameter of around ", diameter_threshold, " mas.")
    list_3_print = input("Print list of those? (y or n) ").lower()
    if list_3_print == "y":
        print(extracted_data3)

    list_3_save = input("Save as .csv? (y or n) ").lower()
    if list_3_save == "y":
        if magnitude == "y":
            output_file = 'stars_visible_bright_big_' + date_str + '.csv'
        else:
            output_file = 'stars_visible_big_' + date_str + '.csv'

        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=extracted_data3.keys())
            writer.writeheader()
            writer.writerows(extracted_data3)

    extracted_data = extracted_data3
