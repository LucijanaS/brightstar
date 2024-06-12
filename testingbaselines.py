"""
Created on: 08.05.2023
Created by: Lucijana Stanic
Function to determine the UV(W) plot for a certain set of telescopes
"""
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

from brightstar_functions import dms_to_decimal, RA_2_HA, R_y, R_x

from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm

import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord


# ---------------------------------------------------------------------
# ---------------------------------------------------------------------


telescopes_VERITAS = [
    {"name": "T1", "latitude": "31 40 29.6N", "longitude": "110 57 03.5W", "height": 1275.893 * u.m},
    {"name": "T2", "latitude": "31 40 28.3N", "longitude": "110 57 06.9W", "height": 1271.016 * u.m},
    {"name": "T3", "latitude": "31 40 32.0N", "longitude": "110 57 07.5W", "height": 1267.358 * u.m},
    {"name": "T4", "latitude": "31 40 30.3N", "longitude": "110 57 10.0W", "height": 1268.273 * u.m},
    {"name": "T5", "latitude": "32 40 30.3N", "longitude": "110 53 10.0W", "height": 1217.273 * u.m},
]

telescopes_VLT = [
    {"name": "T1", "latitude": "24 37 39.4S", "longitude": "70 24 17.7W", "height": 2635 * u.m},
    {"name": "T2", "latitude": "24 37 37.7S", "longitude": "70 24 16.8W", "height": 2635 * u.m},
    {"name": "T3", "latitude": "24 37 36.6S", "longitude": "70 24 15.7W", "height": 2635 * u.m},
    {"name": "T4", "latitude": "24 37 37.2S", "longitude": "70 24 13.8W", "height": 2635 * u.m},
]


telescopes = telescopes_VLT

# UTC offset at the location
utc_offset = -7

# Day (/Night) of observation
date_str = "2019-12-13"


# Extract and convert coordinates
coordinates = []
for telescope in telescopes:
    lat_decimal = dms_to_decimal(telescope['latitude'])
    lon_decimal = dms_to_decimal(telescope['longitude'])
    coordinates.append((lat_decimal, lon_decimal, telescope['height'].value))

# Calculate differences
deltas = {}
for i in range(len(coordinates)):
    for j in range(i + 1, len(coordinates)):
        lat_i, lon_i, height_i = coordinates[i]
        lat_j, lon_j, height_j = coordinates[j]
        delta_lat = lat_j - lat_i
        delta_lon = lon_j - lon_i
        delta_height = height_j - height_i
        deltas[f'delta_lat{j+1}{i+1}'] = delta_lat
        deltas[f'delta_lon{j+1}{i+1}'] = delta_lon
        deltas[f'delta_height{j+1}{i+1}'] = delta_height
        deltas[f'delta_lat{i+1}{j+1}'] = -delta_lat
        deltas[f'delta_lon{i+1}{j+1}'] = -delta_lon
        deltas[f'delta_height{i+1}{j+1}'] = -delta_height


# Earth radius (assuming a spherical Earth)
R = 6371.0 * u.km


# Calculate differences and positional changes
east_diffs = {}
north_diffs = {}
height_diffs = {}

for i in range(len(coordinates)):
    for j in range(i + 1, len(coordinates)):
        lat_i, lon_i, height_i = coordinates[i]
        lat_j, lon_j, height_j = coordinates[j]

        delta_lat = lat_j - lat_i
        delta_lon = lon_j - lon_i
        delta_height = height_j - height_i

        x_E = np.round(((delta_lon * np.pi / 180.0) * R * np.cos(lat_i * np.pi / 180.0)).to(u.meter).value, 3)
        x_N = np.round((delta_lat * np.pi / 180.0 * R).to(u.meter).value, 3)

        east_diffs[f'x_E{j + 1}{i + 1}'] = x_E
        north_diffs[f'x_N{j + 1}{i + 1}'] = x_N
        height_diffs[f'x_up{j + 1}{i + 1}'] = delta_height

        east_diffs[f'x_E{i + 1}{j + 1}'] = -x_E
        north_diffs[f'x_N{i + 1}{j + 1}'] = -x_N
        height_diffs[f'x_up{i + 1}{j + 1}'] = -delta_height

star_of_interest = "β Canis Majoris,Mirzam,0.019,52.632,0.77,1.98,1.75,28000.0,6.378,-16.044,06h 22m 42.0s,-17° 57′ 21″,0.69,0.57,0.52,1.6368e-05,7.967e-06,8.06e-06"

#star_of_interest= "ε Orionis,Alnilam,-0.002,-500.0,0.47,1.7,1.51,30000.0,5.604,-0.798,05h 36m 12.8s,-01° 12′ 07″,0.76,0.63,0.56,2.1577e-05,1.0311e-05,1.0053e-05"
values = star_of_interest.split(',')

BayerF = values[0]
given_ra_decimal = float(values[8])
given_dec_decimal = float(values[9])
diameter_V = float(values[13])
Phi_V = float(values[16])
diameter_in_rad = diameter_V / 1000 * np.pi / (3600 * 180)

# Parse the input string into a datetime object
date_obj = datetime.strptime(date_str, '%Y-%m-%d')

# Convert the date object to a Julian date
date_JD = date_obj.toordinal() + 1721425 + .33333 - (
        1 / 24) * utc_offset  # added 1/3 such since observations will most likely start at 8pm + offset of timezone

# Create a Time object from the observation time in Julian date
observation_time_utc = Time(date_JD, format='jd')

equatorial_coords = SkyCoord(given_ra_decimal, given_dec_decimal, unit=(u.hourangle, u.deg), frame='icrs')

# Define time range for trail calculation
hours_before = -1 / 3600
hours_after = 8.001
start_time = observation_time_utc - TimeDelta(hours_before * u.hour)
end_time = observation_time_utc + TimeDelta(hours_after * u.hour)
times = start_time + (end_time - start_time) * np.linspace(0, 1, 97)[:, None]

datetime_objects = [Time(time[0]).to_datetime() for time in times]

# Extract only the time component from datetime objects and convert to string
time_components = [dt.time().strftime('%H:%M') for dt in datetime_objects]


# Dictionary to hold U, V, W lists
UVW = {}
n_telescopes = len(telescopes)

# Initialize lists for each pair
for i in range(1, n_telescopes+1):
    for j in range(i+1, n_telescopes+1):
        UVW[f'U{i}{j}'] = []
        UVW[f'V{i}{j}'] = []
        UVW[f'W{i}{j}'] = []
        UVW[f'U{j}{i}'] = []
        UVW[f'V{j}{i}'] = []
        UVW[f'W{j}{i}'] = []

# Assumed you have lat1, lat2... defined somewhere, otherwise add this
latitudes = [dms_to_decimal(telescope['latitude']) for telescope in telescopes]

# Iterate over time points
for time in tqdm(times):
    HA_value = RA_2_HA(given_ra_decimal, time.jd)

    # For each telescope pair, calculate UVW
    for i in range(1, n_telescopes + 1):
        lat_i = latitudes[i - 1]  # Adjust for zero-based indexing
        for j in range(i + 1, n_telescopes + 1):
            lat_j = latitudes[j - 1]

            # Retrieve east, north, up differences
            x_E = east_diffs[f'x_E{j}{i}']
            x_N = north_diffs[f'x_N{j}{i}']
            x_up = height_diffs[f'x_up{j}{i}']

            # Apply rotation matrices
            matrices = R_x(given_dec_decimal * u.deg).dot(R_y(HA_value * u.deg)).dot(R_x(-lat_i * u.deg))
            h_plane = np.array([[x_E], [x_N], [x_up]])
            UVW_plane = matrices.dot(h_plane)

            # Append results to UVW dictionary
            UVW[f'U{i}{j}'].append(UVW_plane[0][0])
            UVW[f'V{i}{j}'].append(UVW_plane[1][0])
            UVW[f'W{i}{j}'].append(UVW_plane[2][0])
            # Append the reverse as well
            UVW[f'U{j}{i}'].append(-UVW_plane[0][0])
            UVW[f'V{j}{i}'].append(-UVW_plane[1][0])
            UVW[f'W{j}{i}'].append(-UVW_plane[2][0])

# Calculate Mlambda for the wavelength
wavelength = 540e-9  # wavelength in meters
Mlambda = wavelength * 1e6  # convert meters to microns

# Track labels to prevent duplicates in the legend
plotted_labels = set()

plt.figure(figsize=(10, 8))
for i in range(1, n_telescopes + 1):
    for j in range(i + 1, n_telescopes + 1):
        label = f"T{i}T{j}"
        plt.plot(np.array(UVW[f'U{i}{j}']) / Mlambda, np.array(UVW[f'V{i}{j}']) / Mlambda, '.', label=label)
        plt.plot(np.array(UVW[f'U{j}{i}']) / Mlambda, np.array(UVW[f'V{j}{i}']) / Mlambda, '.',
                     color=plt.gca().lines[-1].get_color())  # Match color with the previous plot
        plotted_labels.add(label)

axis = 300
plt.xlim(-axis, axis)
plt.ylim(-axis, axis)
plt.xlabel("u ($M\lambda$)")
plt.ylabel("v ($M\lambda$)")
plt.legend()
plt.title(f"{BayerF} {date_str}")
plt.show()