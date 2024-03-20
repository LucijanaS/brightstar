"""
Created on: 07.03.2023
Created by: Lucijana Stanic

"""
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import astropy.units as u
from astropy.time import Time, TimeDelta
from datetime import date, datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from scipy.special import j1, j0



# - * - coding: utf - 8 - * -

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

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


date_str = "2024-05-22"

utc_offset = 0

lat_deg = "28 17 58.8N"  # Telescopio Carlos Sanchez
lon_deg = "16 30 39.7E"

lat = dms_to_decimal(lat_deg)
lon = dms_to_decimal(lon_deg)

height = 2381.25 * u.m


# Parse the input string into a datetime object
date_obj = datetime.strptime(date_str, '%Y-%m-%d')

# Convert the date object to a Julian date
utc_offset = 1
date_JD = date_obj.toordinal() + 1721425 + .33333 - (
        1 / 24) * utc_offset  # added 1/3 such since observations will most likely start at 8pm + offset of timezone

# Create a Time object from the observation time in Julian date
observation_time_utc = Time(date_JD, format='jd')

# Given RA and Dec
given_ra = "13h 25m 11.6s"
given_dec = "-11° 09′ 41″"


# Convert given RA and Dec to the format in the JSON data for comparison
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


given_ra_decimal, given_dec_decimal = convert_ra_dec(given_ra, given_dec)


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
        AltAz(obstime=time, location=EarthLocation(lat=lat, lon=lon, height=height)))
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
plt.title('Star Trail')
plt.xlabel('Time')
plt.ylabel('Altitude [°]')
plt.ylim(0, 90)
plt.grid(True)
plt.show()


def R_x(a):
    return np.array([[1, 0, 0],
                     [0, np.cos(a), -np.sin(a)],
                     [0, np.sin(a), np.cos(a)]])


def R_y(b):
    return np.array([[np.cos(b), 0, np.sin(b)],
                     [0, 1, 0],
                     [-np.sin(b), 0, np.cos(b)]])


def RA_2_HA(right_ascension, local_time):
    """Converts right ascension (in decimal degrees) to hour angle (in decimal degrees) given the observation time (Julian date)."""
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


U = []
V = []
W = []


# Create a grid of points
resolution = 300
size_to_plot = 100
x = np.linspace(-size_to_plot, size_to_plot, resolution)
y = np.linspace(-size_to_plot, size_to_plot, resolution)
X, Y = np.meshgrid(x, y)
R = np.sqrt(X**2 + Y**2)


for time in times:
    HA_value = RA_2_HA(given_ra_decimal, time.jd)
    matrices = R_x(given_dec_decimal * u.deg).dot(R_y(HA_value * u.deg)).dot(R_x(-lat * u.deg))
    diff_N = -92.5
    diff_E = -13
    diff_H = 5.5
    h_plane = np.array([[diff_N], [diff_E], [diff_H]])
    UVW_plane = matrices.dot(h_plane)

    U.append(UVW_plane[0][0])  # Append the first element of the first dimension
    V.append(UVW_plane[1][0])  # Append the first element of the second dimension
    W.append(UVW_plane[2][0])  # Append the first element of the third dimension

resolution = 300
diameter = 0.01
diameter_in_rad = diameter*np.pi/(3600*180)
diameter_mas = diameter*1000
wavelength = (5.4e-7) # wavelength in meters
wavelength_nm = wavelength *10**9
A = calculate_covered_area(U, V)
intensity_values = intensity(R, diameter_in_rad, wavelength)
plt.imshow(intensity_values, extent=(-size_to_plot, size_to_plot, -size_to_plot, size_to_plot), origin='lower')
plt.plot(U, V, '.')
plt.gca().set_aspect('equal')
plt.show()
