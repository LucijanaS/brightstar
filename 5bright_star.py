"""
Created on: 07.03.2023
Created by: Lucijana Stanic

By entering a row from one of the bright star files, this script returns two plots for that particular star.
One plot is its path across the night sky during the chosen night while one shows the UV-path it traces
"""
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

from brightstar_input import lat_deg1, lon_deg1, height1, date_str, x_up, x_E, x_N, utc_offset
from brightstar_functions import dms_to_decimal, convert_ra_dec, RA_2_HA, R_y, R_x, calculate_covered_area, visibility

from datetime import datetime, timedelta

import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

# - * - coding: utf - 8 - * -

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
"""
# Given RA and Dec
given_ra = "16h 27m 01.4s"
given_dec = "-18° 27′ 23″"
given_ra_decimal, given_dec_decimal = convert_ra_dec(given_ra, given_dec)

# Given diameter in milliarcsecond
diameter_mas = 0.19
diameter_in_rad = diameter_mas / 1000 * np.pi / (3600 * 180)
"""
star_of_interest = input("Star of interest (input row from .csv file):")
#star_of_interest = "α Boötis,Arcturus,0.09,11.111,2.46,-0.04,1.19,4850.0,14.261,19.183,14h 15m 39.7s,+19° 10′ 57″,9.26,14.02,10.84,3.451e-06,5.1202e-05,1.3499e-05,4456.530301466944,0.0016"
values = star_of_interest.split(',')

BayerF = values[0]
given_ra_decimal = float(values[8])
given_dec_decimal = float(values[9])
diameter_V = float(values[13])
Phi_V = float(values[16])
diameter_in_rad = diameter_V / 1000 * np.pi / (3600 * 180)

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

lat = dms_to_decimal(lat_deg1)
lon = dms_to_decimal(lon_deg1)

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
hours_after = 12.001
start_time = observation_time_utc - TimeDelta(hours_before * u.hour)
end_time = observation_time_utc + TimeDelta(hours_after * u.hour)
times = start_time + (end_time - start_time) * np.linspace(0, 1, 97)[:, None]

# Calculate coordinates at each time step
altitudes = []
azimuths = []

for time in times:
    altaz_coords = equatorial_coords.transform_to(
        AltAz(obstime=time, location=EarthLocation(lat=lat, lon=lon, height=height1)))
    altitude = altaz_coords.alt
    azimuth = altaz_coords.az
    altitudes.append(altitude)
    azimuths.append(azimuth)

# Convert lists to arrays
altitudes = np.array(altitudes)
azimuths = np.array(azimuths)
azimuths_flat = azimuths.flatten()
datetime_objects = [Time(time[0]).to_datetime() for time in times]

# Extract only the time component from datetime objects and convert to string
time_components = [dt.time().strftime('%H:%M') for dt in datetime_objects]

# Plot the trail
plt.scatter(time_components, altitudes, c=azimuths_flat)
plt.colorbar(label='Azimuth [°]')  # Add color bar indicating azimuth
plt.xticks(time_components[::16], rotation=0)
plt.title(BayerF)
plt.xlabel('Time')
plt.ylabel('Altitude [°]')
plt.ylim(0, 90)
plt.grid(True)
plt.show()
plt.clf()

U = []
V = []
W = []

# Create a grid of points
resolution = 300
size_to_plot = np.sqrt(x_E**2+x_N**2)
x = np.linspace(-size_to_plot, size_to_plot, resolution)
y = np.linspace(-size_to_plot, size_to_plot, resolution)
X, Y = np.meshgrid(x, y)
R = np.sqrt(X ** 2 + Y ** 2)

for time in times:
    HA_value = RA_2_HA(given_ra_decimal, time.jd)
    matrices = R_x(given_dec_decimal * u.deg).dot(R_y(HA_value * u.deg)).dot(R_x(-lat * u.deg))
    h_plane = np.array([[x_E], [x_N], [x_up]])
    UVW_plane = matrices.dot(h_plane)

    U.append(UVW_plane[0][0])
    V.append(UVW_plane[1][0])
    W.append(UVW_plane[2][0])

resolution = 300
wavelength = (5.4e-7)  # wavelength in meters
wavelength_nm = wavelength * 10 ** 9
A = calculate_covered_area(U, V)
intensity_values = visibility(R, diameter_in_rad, wavelength)
plt.imshow(intensity_values, norm=None, extent=(-size_to_plot, size_to_plot, -size_to_plot, size_to_plot), origin='lower',
           cmap='gray')
plt.plot(U, V, '.', color='gold', markeredgecolor='black')
plt.title(BayerF + " diameter: " + str(diameter_V) + " mas\n "
                                                     "$\Phi$ = " + str(
    np.round(Phi_V, 7)) + " photons m$^{-2}$ s$^{-1}$ Hz$^{-1}$")
plt.xlabel('U [m]')
plt.ylabel('V [m]')
plt.colorbar(label='Intensity')
plt.gca().set_aspect('equal')
plt.show()
plt.clf()
