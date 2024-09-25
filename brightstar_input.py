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


# Location

# Main observation coordinates, input: "dd mm ss D"
#lat_deg1 = "43 45 03.5N"   # Crete SKINAKAS 1.3
#lon_deg1 = "6 55 16.2E"
#height1 = 1750 * u.m

# Second observation point to determine path on UVW plane
#lat_deg2 = "43 45 02.9N"  # Crete SKINAKAS 0.6
#lon_deg2 = "6 55 17.2E"
#height2 = 1750 * u.m

lat_deg1 = "23 20 31.9N"   # gamsburg
lon_deg1 = "16 13 29.7E"
height1 = 1750 * u.m

"""
# Main observation coordinates, input: "dd mm ss D"
lat_deg1 = "28 18 01.8N"  # Telescopio Carlos Sánchez
lon_deg1 = "16 30 39.2W"
height1 = 2386.75 * u.m

# Second observation point to determine path on UVW plane
lat_deg2 = "28 17 58.8N"  # IAC80 Telescope
lon_deg2 = "16 30 39.7W"
height2 = 2381.25 * u.m


lat_deg1 = "47 23 51.9N"  # Irchel Zurich
lon_deg1 = "8 32 53.9E"
height1 = 400 * u.m


lat_deg1 = "46 13 42N"  # St. Luc
lon_deg1 = "7 36 45E"
height1 = 2176 * u.m
"""

# Convert coordinates to decimal degrees
#lat1 = dms_to_decimal(lat_deg1)
#lon1 = dms_to_decimal(lon_deg1)
#lat2 = dms_to_decimal(lat_deg2)
#lon2 = dms_to_decimal(lon_deg2)

lat1 = 43.75370
lon1 = 6.92294
lat2 = 43.75370
lon2 = 6.92312


# Calculate differences
delta_lat = lat2 - lat1
delta_lon = lon2 - lon1
delta_height = 0.0

# Earth radius (assuming a spherical Earth)
R = 6371.0 * u.km





# Calculate differences in east (x_E) and north (x_N) directions
x_E = np.round(((delta_lon * np.pi / 180.0) * R * np.cos(lat1 * np.pi / 180.0)).value*1000, 3)
x_N = np.round((delta_lat * R * np.pi / 180.0).value*1000, 3)

# Difference in height (x_up)
#x_up = delta_height.value
x_up = 0

"""
x_E = 119
x_N = 120
x_up = 0
"""
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

# Time

# UTC offset at the location
utc_offset = +2

# Day (/Night) of observation
date_str = "2024-10-15"

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

# How many stars of the brightest stars should be extracted and taken into account?
n_brightest_stars = 100
#print(np.sqrt(x_E**2 + x_N**2))