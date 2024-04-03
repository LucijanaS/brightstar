"""
Created on: 13.03.2023
Created by: Lucijana Stanic

This script calculates the UVW-trace for a certain location and time input in brightstar_input.py and determines values
from that that might be useful in determining which stars could be the best candidates for observations.
"""
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

from brightstar_functions import dms_to_decimal, RA_2_HA, R_y, R_x, calculate_covered_area, intensity
from brightstar_input import lat_deg1, utc_offset, date_str, x_E, x_N, x_up

import csv
from datetime import datetime
from tqdm import tqdm

import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord

# - * - coding: utf - 8 - * -

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

# Initialize a dictionary to store data from each column
data = {}
file_name = "stars_visible_"+date_str+".csv"
# Open the CSV file
with open(file_name, newline='') as csvfile:
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

# Convert coordinates from deg to dec
lat_dec1 = dms_to_decimal(lat_deg1)

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

# Evaluate motion on the UVW plane of each star on a night
RA = []
for star in data["RA_decimal"]:
    RA.append(star)

Dec = []
for star in data["Dec_decimal"]:
    Dec.append(star)

diameter_V = []
for star in data["Diameter_V"]:
    diameter_V.append(star)

equatorial_coords = []
for i in range(len(data["RA_decimal"])):
    coord = SkyCoord(RA[i], Dec[i], unit=(u.hourangle, u.deg), frame='icrs')
    equatorial_coords.append(coord)

# Initialize lists to store altitude and azimuth values for each star
U_per_star = []
V_per_star = []
W_per_star = []
intensity_sum_per_star = []
intensity_map_per_star = []
intensity_std_per_star = []
intensity_max_min_per_star = []
covered_area = []
# Iterate over each equatorial coordinate
n = len(equatorial_coords)

# Create a grid of points
resolution = 300
size_to_plot = 30
x = np.linspace(-size_to_plot, size_to_plot, resolution)
y = np.linspace(-size_to_plot, size_to_plot, resolution)
X, Y = np.meshgrid(x, y)
R = np.sqrt(X ** 2 + Y ** 2)

print("Determining intensity variations and covered area in UVW-plane for each star in", file_name, "during the night of", date_str)
for i in tqdm(range(0, n)):
    # Get RA and Dec for the current star using its index i
    given_ra_decimal = float(RA[i])
    given_dec_decimal = float(Dec[i])
    current_diameter_V = float(diameter_V[i])
    diameter_in_rad = current_diameter_V / 1000 * np.pi / (3600 * 180)
    wavelength = 540 * 10 ** -9
    U = []
    V = []
    W = []
    intensity_value = []

    # Define the time span
    hours_before = 0
    hours_after = 12
    start_time = observation_time_utc - TimeDelta(hours_before * u.hour)
    end_time = observation_time_utc + TimeDelta(hours_after * u.hour)
    intensity_map = intensity(R, diameter_in_rad, wavelength)

    # Calculate altitude and azimuth for each time point
    times = start_time + (end_time - start_time) * np.linspace(0, 1, 97)[:, None]
    for time in times:
        HA_value = RA_2_HA(given_ra_decimal, time.jd)
        matrices = R_x(given_dec_decimal * u.deg).dot(R_y(HA_value * u.deg)).dot(R_x(-lat_dec1 * u.deg))
        h_plane = np.array([[x_N], [x_E], [x_up]])
        UVW_plane = matrices.dot(h_plane)

        U.append(UVW_plane[0][0])
        V.append(UVW_plane[1][0])
        W.append(UVW_plane[2][0])
        intensity_value.append(intensity_map[int(UVW_plane[0][0]), int(UVW_plane[1][0])])

    plt.plot(U, V, '.')
    plt.gca().set_aspect('equal')
    plt.title(data['BayerF'][i])
    plt.imshow(intensity_map, extent=(-size_to_plot, size_to_plot, -size_to_plot, size_to_plot), origin='lower')
    plt.colorbar(label='Intensity')

    # plt.show()
    plt.clf()
    A = calculate_covered_area(U, V)

    U_per_star.append(U)
    V_per_star.append(V)
    W_per_star.append(W)
    intensity_map_per_star.append(intensity_map)
    intensity_max_min_per_star.append(np.max(intensity_value)-np.min(intensity_value))
    intensity_sum_per_star.append(np.round(np.sum(intensity_value), 4))
    intensity_std_per_star.append(np.std(intensity_value))
    covered_area.append(A)

# Add the new column of data to the dictionary
data["Area"] = covered_area
data["Intensity Sum"] = intensity_sum_per_star
data["Intensity Std"] = intensity_std_per_star
data["Intensity Max-Min"] = intensity_max_min_per_star

list_save = "y"
if list_save == "y":
    output_file = file_name

    # Write the extracted data to the CSV file
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Write the header row using the keys of the extracted_data dictionary
        writer.writerow(data.keys())

        # Write the data rows
        # Determine the number of rows based on the length of one of the lists in extracted_data
        num_rows = len(next(iter(data.values())))

        for i in range(num_rows):
            # Get the data for each column and write it to the CSV file
            row_data = [data[key][i] for key in data.keys()]
            writer.writerow(row_data)
    print("Added to", output_file)

indices_sorted_high = sorted(range(len(data['Intensity Sum'])), key=lambda i: data['Intensity Max-Min'][i],
                             reverse=True)[:5]


#top_stars = input("Show plots of the stars with highest intensity variations? (y/n)")
top_stars = "y"
if top_stars == "y":
    for i in indices_sorted_high:
        plt.plot(U_per_star[i], V_per_star[i], '.', color='gold', markeredgecolor='black')
        plt.gca().set_aspect('equal')
        plt.title(data['BayerF'][i], )
        plt.imshow(intensity_map_per_star[i], cmap='gray',
               extent=(-size_to_plot, size_to_plot, -size_to_plot, size_to_plot), origin='lower')
        plt.colorbar(label='Intensity')

        plt.show()
        plt.clf()


indices_sorted_low = sorted(range(len(data['Intensity Sum'])), key=lambda i: data['Intensity Max-Min'][i])[:5]

for i in indices_sorted_low:
    plt.plot(U_per_star[i], V_per_star[i], '.', color='gold', markeredgecolor='black')
    plt.gca().set_aspect('equal')
    plt.title(data['BayerF'][i], )
    plt.imshow(intensity_map_per_star[i], cmap='gray',
               extent=(-size_to_plot, size_to_plot, -size_to_plot, size_to_plot), origin='lower')
    plt.colorbar(label='Intensity')

    plt.show()
    plt.clf()
