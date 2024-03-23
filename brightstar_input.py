"""
Created on: 20.03.2023
Created by: Lucijana Stanic

Here one inputs the date and location of the observation, this is input is then used in all the scripts
"""
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

from brightstar_functions import dms_to_decimal

import astropy.units as u
from astropy.coordinates import Angle
import numpy as np


# ---------------------------------------------------------------------
# ---------------------------------------------------------------------


# Main observation coordinates, input: "dd mm ss D"
lat_deg1 = "28 18 01.8N"  # Telescopio Carlos Sánchez
lon_deg1 = "16 30 39.2W"
height1 = 2386.75 * u.m

# Second observation point to determine path on UVW plane
lat_deg2 = "28 17 58.8N"  # IAC80 Telescope
lon_deg2 = "16 30 39.7W"
height2 = 2381.25 * u.m

# UTC offset at the location
utc_offset = 0

# Day (/Night) of observation
date_str = "2024-05-22"

# How many stars of the brightest stars should be extracted and taken into account?
n_brightest_stars = 1000

# Convert coordinates to decimal degrees
lat1 = dms_to_decimal(lat_deg1)
lon1 = dms_to_decimal(lon_deg1)
lat2 = dms_to_decimal(lat_deg2)
lon2 = dms_to_decimal(lon_deg2)

# Calculate differences
delta_lat = lat2 - lat1
delta_lon = lon2 - lon1
delta_height = height2 - height1

# Earth radius (assuming a spherical Earth)
R = 6371.0 * u.km

# Calculate differences in east (x_E) and north (x_N) directions
x_E = np.round(((delta_lon * np.pi / 180.0) * R * np.cos(lat1 * np.pi / 180.0)).value*1000, 3)
x_N = np.round((delta_lat * R * np.pi / 180.0).value*1000, 3)

# Difference in height (x_up)
x_up = delta_height.value