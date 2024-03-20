"""
Created on: 13.03.2023
Created by: Lucijana Stanic

"""
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

import numpy as np
import astropy.units as u
import csv
from datetime import datetime
from astropy.coordinates import SkyCoord
from astropy.time import Time, TimeDelta
from tqdm import tqdm
import matplotlib.pyplot as plt
import scipy.signal as signal
from scipy.special import j1, j0
from scipy.spatial import ConvexHull

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


def dms_to_decimal(dms_str):
    # Split the string into degrees, minutes, seconds, and direction
    degrees, minutes, seconds = map(float, dms_str[:-1].split(' '))

    # Extract direction
    direction = dms_str[-1]

    # Calculate the decimal degrees
    decimal_degrees = degrees + (minutes / 60.0) + (seconds / 3600.0)

    # Adjust for negative values if direction is South or West
    if direction in ['S', 'W']:
        decimal_degrees *= -1

    return decimal_degrees


def R_x(a):
    return np.array([[1, 0, 0],
                     [0, np.cos(a), -np.sin(a)],
                     [0, np.sin(a), np.cos(a)]])


def R_y(b):
    return np.array([[np.cos(b), 0, np.sin(b)],
                     [0, 1, 0],
                     [-np.sin(b), 0, np.cos(b)]])


def RA_2_HA(right_ascension, local_time):
    """Converts right ascension (in decimal degrees) to hour angle (in decimal degrees) given the observation time (
    Julian date)."""
    # Calculate Greenwich Mean Sidereal Time (GMST) at 0h UTC on the given observation date
    # GMST at 0h on January 1, 2000 (J2000 epoch) is 280.4606°
    gmst_0h_J2000 = 280.4606

    # Calculate the number of Julian centuries since J2000 epoch
    T = (local_time - 2451545.0) / 36525.0

    # Calculate the GMST at the given observation time
    gmst = gmst_0h_J2000 + 360.98564724 * (local_time - 2451545.0) + 0.000387933 * T ** 2 - (T ** 3) / 38710000.0

    # Normalize GMST to the range [0, 360) degrees
    gmst %= 360.0

    # Calculate the hour angle (HA) in decimal degrees
    ha = gmst - right_ascension

    # Normalize hour angle to the range [-180, 180) degrees
    ha = (ha + 180.0) % 360.0 - 180.0

    return float(ha[0])


def convert_ra_dec(ra_str, dec_str):
    ra_parts = ra_str.split(' ')
    ra_h = int(ra_parts[0][:-1])
    ra_m = int(ra_parts[1][:-1])
    ra_s = float(ra_parts[2][:-1])
    ra_decimal = ra_h + ra_m / 60 + ra_s / 3600

    dec_parts = dec_str.split(' ')
    dec_d = int(dec_parts[0][:-1])
    dec_m = int(dec_parts[1][:-1])
    dec_s = float(dec_parts[2][:-1])
    dec_decimal = dec_d + dec_m / 60 + dec_s / 3600

    return ra_decimal, dec_decimal


def calculate_covered_area(U, V):
    # Combine U and V coordinates into a single array
    points = np.column_stack((U, V))

    # Calculate the convex hull of the points
    hull = ConvexHull(points)

    # Calculate the area of the convex hull
    area = hull.volume if len(U) > 2 else 0

    return area

def intensity(b, theta, lambda_):
    input = np.pi * b * theta / lambda_
    B_1 = j1(input)
    return B_1


# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

# Input date and location (either in degrees or decimal values, choose which input to use)

date_str = "2024-05-22"

utc_offset = 0

height = 2381.25 * u.m

lat_deg = "28 17 58.8N"  # Telescopio Carlos Sanchez
lon_deg = "16 30 39.7E"

lat_deg_to_dec = dms_to_decimal(lat_deg)
lon_deg_to_dec = dms_to_decimal(lon_deg)

lat_dec = 47.3769 * u.deg
lon_dec = 8.5417 * u.deg

loc_input = 'deg'

if loc_input == 'dec':
    lat = lat_dec
    lon = lon_dec
if loc_input == 'deg':
    lat = lat_deg_to_dec
    lon = lon_deg_to_dec

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






for i in tqdm(range(0, n)):
    # Get RA and Dec for the current star using its index i
    given_ra_decimal = float(RA[i])
    given_dec_decimal = float(Dec[i])
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
        matrices = R_x(given_dec_decimal * u.deg).dot(R_y(HA_value * u.deg)).dot(R_x(-lat * u.deg))
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
    intensity_values = intensity(R, diameter, wavelength)

    # plt.show()
    A = calculate_covered_area(U, V)

    U_per_star.append(U)
    V_per_star.append(V)
    W_per_star.append(W)
    covered_area.append(A)




"""

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

print(max(covered_area))
"""