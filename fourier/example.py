import numpy as np
import matplotlib.pyplot as plt
import sys

from fourier import grids
from sources import blob, smooth, crescent, eclipsed_sphere, star_with_disk
from visib import sbright, correldens
from graphics import draw

import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from datetime import datetime

import matplotlib.colors as mcolors


def parse_colormap(file_path):
    temperatures = []
    colors = []
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            temperature = int(parts[0])
            color = parts[2]
            temperatures.append(temperature)
            colors.append(color)
    return temperatures, colors

def create_custom_colormap(temperatures, colors):
    norm = mcolors.Normalize(vmin=min(temperatures), vmax=max(temperatures))
    tuples = list(zip(map(norm, temperatures), colors))
    cmap = mcolors.LinearSegmentedColormap.from_list("custom_cmap", tuples)
    return cmap


temperatures, bb_colors = parse_colormap('blackbody_colors')
bb_cmap = create_custom_colormap(temperatures, bb_colors)


plt.style.use('dark_background')

N = 2048  # N is grid size (equal on ground and sky)
ds = 1e-10  # ds is grid spacing on the sky

lam = 500e-9 # wavelength in meters

sx, sy, x, y = grids(ds, N, lam)

Teff = 1e4
temp = 27 # for sensible UV image
temp = 0.6 # for correct phi on the plot
temp = 0.35 # for correct sum of phi
temp = 1.3 # for correct mag

tempratio = 0.00001
Tmap = Teff * star_with_disk(sx, sy, core_radius=0.22e-8, disk_radius=2e-8, elong=1.49, core_temp=0.4, disk_temp=0.13, pa=(20/180)*np.pi)
Tmap = smooth(Tmap, 3)
draw(sx, sy, Tmap, 1, 'sky', cmap=bb_cmap, title='$T_{\\rm eff}$')
plt.show()

S = sbright(Tmap, lam, ds)
f = correldens(S, lam)

#draw(x, y, f, 50, 'ground', cmap='gray', title=r'$\Phi\,|V|^2$')
#plt.show()


def parse_colormap(file_path):
    temperatures = []
    colors = []
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            temperature = int(parts[0])
            color = parts[2]
            temperatures.append(temperature)
            colors.append(color)
    return temperatures, colors

def create_custom_colormap(temperatures, colors):
    norm = mcolors.Normalize(vmin=min(temperatures), vmax=max(temperatures))
    tuples = list(zip(map(norm, temperatures), colors))
    cmap = mcolors.LinearSegmentedColormap.from_list("custom_cmap", tuples)
    return cmap



#star_of_interest = input("Star of interest (input row from .csv file):")
star_of_interest = "γ Cassiopeiae,Navi,0.016,62.5,1.24,2.47,2.32,30000.0,0.945,60.717,00h 56m 42.5s,+60° 43′ 00″,0.53,0.44,0.39,1.0616824921819226e-05,5.073258968632882e-06,4.767772686428327e-06"
values = star_of_interest.split(',')



BayerF = values[0]
given_ra_decimal = float(values[8])
given_dec_decimal = float(values[9])
diameter_V = float(values[13])
Phi_V = float(values[16])
diameter_in_rad = diameter_V / 1000 * np.pi / (3600 * 180)

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

lat1 = 43.75370
lon1 = 6.92294
lat2 = 43.75370
lon2 = 6.92312


# Calculate differences
delta_lat = lat2 - lat1
delta_lon = lon2 - lon1
height1=0.0
delta_height = 0.0
# Earth radius (assuming a spherical Earth)
R = 6371.0 * u.km
x_E = np.round(((delta_lon * np.pi / 180.0) * R * np.cos(lat1 * np.pi / 180.0)).value*1000, 3)
x_N = np.round((delta_lat * R * np.pi / 180.0).value*1000, 3)
x_up = 0


# Earth radius (assuming a spherical Earth)
R = 6371.0 * u.km

utc_offset = 1


date_str = "2024-11-04"

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


U = []
V = []
W = []



def R_x(a):
    """Part of a rotation matrix, split into R_x and R_y for better overview of what is happening"""
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

for time in times:
    HA_value = RA_2_HA(given_ra_decimal, time.jd)
    matrices = R_x(given_dec_decimal * u.deg).dot(R_y(HA_value * u.deg)).dot(R_x(-lat1 * u.deg))
    h_plane = np.array([[x_E], [x_N], [x_up]])
    UVW_plane = matrices.dot(h_plane)

    U.append(UVW_plane[0][0])
    V.append(UVW_plane[1][0])
    W.append(UVW_plane[2][0])

resolution = 300
wavelength = (5.4e-7)  # wavelength in meters
wavelength_nm = wavelength * 10 ** 9
plt.plot(U, V, '.', color='gold', markeredgecolor='black')
draw(x, y, f, 50, 'ground', cmap='gray', title='$\Phi\,|V|^2$')
plt.plot(U, V, '.', color='gold', markeredgecolor='black')
U_neg = [-u for u in U]
V_neg = [-v for v in V]

plt.plot(U_neg, V_neg, '.', color='gold', markeredgecolor='black')


plt.xlabel('U [m]')
plt.ylabel('V [m]')
#plt.colorbar(label='Intensity')
plt.gca().set_aspect('equal')
plt.show()
plt.clf()

print(np.sum(f))