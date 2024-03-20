"""
Created on: 13.03.2023
Created by: Lucijana Stanic

"""
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

from brightstar_functions import dms_to_decimal, RA_2_HA, R_y, R_x, calculate_covered_area
from brightstar_input import lat_deg1, utc_offset, date_str

import numpy as np
import astropy.units as u
import csv
from datetime import datetime
from astropy.coordinates import SkyCoord
from astropy.time import Time, TimeDelta
from tqdm import tqdm
import matplotlib.pyplot as plt
from scipy.special import j1


# - * - coding: utf - 8 - * -

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

# Initialize a dictionary to store data from each column
data = {}

# Open the CSV file
with open('stars_visible_2024-05-22.csv', newline='') as csvfile:
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

def intensity(b, theta, lambda_):
    input = np.pi * b * theta / lambda_
    B_1 = j1(input)
    return B_1

# Input date and location (either in degrees or decimal values, choose which input to use)
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
covered_area = []
# Iterate over each equatorial coordinate
n = len(equatorial_coords)



# Create a grid of points
resolution = 300
size_to_plot = 100
x = np.linspace(-size_to_plot, size_to_plot, resolution)
y = np.linspace(-size_to_plot, size_to_plot, resolution)
X, Y = np.meshgrid(x, y)
R = np.sqrt(X**2 + Y**2)


for i in tqdm(range(0, n)):
    # Get RA and Dec for the current star using its index i
    given_ra_decimal = float(RA[i])
    given_dec_decimal = float(Dec[i])
    current_diameter_V = float(diameter_V[i][0])
    diameter_in_rad = current_diameter_V/1000 * np.pi / (3600 * 180)
    wavelength = 540*10**-9
    U = []
    V = []
    W = []

    # Define the time span
    hours_before = 0
    hours_after = 12
    start_time = observation_time_utc - TimeDelta(hours_before * u.hour)
    end_time = observation_time_utc + TimeDelta(hours_after * u.hour)

    # Calculate altitude and azimuth for each time point
    times = start_time + (end_time - start_time) * np.linspace(0, 1, 97)[:, None]
    for time in times:
        HA_value = RA_2_HA(given_ra_decimal, time.jd)
        matrices = R_x(given_dec_decimal * u.deg).dot(R_y(HA_value * u.deg)).dot(R_x(-lat_dec1 * u.deg))
        diff_N = -92.5
        diff_E = -13
        diff_H = 5.5
        h_plane = np.array([[diff_N], [diff_E], [diff_H]])
        UVW_plane = matrices.dot(h_plane)

        U.append(UVW_plane[0][0])
        V.append(UVW_plane[1][0])
        W.append(UVW_plane[2][0])

    plt.plot(U, V, '.')
    plt.gca().set_aspect('equal')
    plt.title(data['BayerF'][i])
    intensity_values = intensity(R, diameter_in_rad, wavelength)
    #plt.imshow(intensity_values, extent=(-size_to_plot, size_to_plot, -size_to_plot, size_to_plot), origin='lower')

    #plt.show()
    A = calculate_covered_area(U, V)

    U_per_star.append(U)
    V_per_star.append(V)
    W_per_star.append(W)
    covered_area.append(A)





# Add the new column of data to the dictionary
data["Area"] = covered_area
list_save = input("Save as .csv? (y or n) ")
if list_save == "y":
    output_file = 'stars_visible_' + date_str + '_area.csv'

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
    print("Saved as", output_file)
