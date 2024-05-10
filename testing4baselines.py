"""
Created on: 27.03.2023
Created by: Lucijana Stanic
Python file created to test the UVW coverage determined and cross check it with VERTIAS' results
"""
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

from brightstar_functions import dms_to_decimal, convert_ra_dec, RA_2_HA, R_y, R_x, calculate_covered_area, visibility

from datetime import datetime, timedelta

import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm

import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, EarthLocation, AltAz


# ---------------------------------------------------------------------
# ---------------------------------------------------------------------


# Main observation coordinates, input: "dd mm ss D"
lat_deg1 = "31 40 29.6N"
lon_deg1 = "110 57 03.5W"
height1 = 1275.893 * u.m

# Second observation point to determine path on UVW plane
lat_deg2 = "31 40 28.3N"
lon_deg2 = "110 57 06.9W"
height2 = 1271.016 * u.m

# Second observation point to determine path on UVW plane
lat_deg3 = "31 40 32.0N"
lon_deg3 = "110 57 07.5W"
height3 = 1267.358 * u.m

# Second observation point to determine path on UVW plane
lat_deg4 = "31 40 30.3N"
lon_deg4 = "110 57 10.0W"
height4 = 1268.273 * u.m

# UTC offset at the location
utc_offset = -7

# Day (/Night) of observation
date_str = "2019-12-13"


# Convert coordinates to decimal degrees
lat1 = dms_to_decimal(lat_deg1)
lon1 = dms_to_decimal(lon_deg1)
lat2 = dms_to_decimal(lat_deg2)
lon2 = dms_to_decimal(lon_deg2)
lat3 = dms_to_decimal(lat_deg3)
lon3 = dms_to_decimal(lon_deg3)
lat4 = dms_to_decimal(lat_deg4)
lon4 = dms_to_decimal(lon_deg4)

# Calculate differences
delta_lat21 = lat2 - lat1
delta_lat12 = lat1 - lat2

delta_lon21 = lon2 - lon1
delta_lon12 = lon1 - lon2

delta_height21 = height2 - height1
delta_height12 = height1 - height2


delta_lat31 = lat3 - lat1
delta_lat13 = lat1 - lat3

delta_lon31 = lon3 - lon1
delta_lon13 = lon1 - lon3

delta_height31 = height3 - height1
delta_height13 = height1 - height3


delta_lat41 = lat4 - lat1
delta_lat14 = lat1 - lat4

delta_lon41 = lon4 - lon1
delta_lon14 = lon1 - lon4

delta_height41 = height4 - height1
delta_height14 = height1 - height4


delta_lat32 = lat3 - lat2
delta_lat23 = lat2 - lat3

delta_lon32 = lon3 - lon2
delta_lon23 = lon2 - lon3

delta_height32 = height3 - height2
delta_height23 = height2 - height3


delta_lat42 = lat4 - lat2
delta_lat24 = lat2 - lat4

delta_lon42 = lon4 - lon2
delta_lon24 = lon2 - lon4

delta_height42 = height4 - height2
delta_height24 = height2 - height4


delta_lat43 = lat4 - lat3
delta_lat34 = lat3 - lat4

delta_lon43 = lon4 - lon3
delta_lon34 = lon3 - lon4

delta_height43 = height4 - height3
delta_height34 = height3 - height4


# Earth radius (assuming a spherical Earth)
R = 6371.0 * u.km

# Calculate differences in east (x_E) and north (x_N) directions
x_E12 = np.round(((delta_lon12 * np.pi / 180.0) * R * np.cos(lat1 * np.pi / 180.0)).value*1000, 3)
x_E21 = np.round(((delta_lon21 * np.pi / 180.0) * R * np.cos(lat1 * np.pi / 180.0)).value*1000, 3)
x_E13 = np.round(((delta_lon13 * np.pi / 180.0) * R * np.cos(lat1 * np.pi / 180.0)).value*1000, 3)
x_E31 = np.round(((delta_lon31 * np.pi / 180.0) * R * np.cos(lat1 * np.pi / 180.0)).value*1000, 3)
x_E14 = np.round(((delta_lon14 * np.pi / 180.0) * R * np.cos(lat1 * np.pi / 180.0)).value*1000, 3)
x_E41 = np.round(((delta_lon41 * np.pi / 180.0) * R * np.cos(lat1 * np.pi / 180.0)).value*1000, 3)
x_E23 = np.round(((delta_lon23 * np.pi / 180.0) * R * np.cos(lat1 * np.pi / 180.0)).value*1000, 3)
x_E32 = np.round(((delta_lon32 * np.pi / 180.0) * R * np.cos(lat1 * np.pi / 180.0)).value*1000, 3)
x_E24 = np.round(((delta_lon24 * np.pi / 180.0) * R * np.cos(lat1 * np.pi / 180.0)).value*1000, 3)
x_E42 = np.round(((delta_lon42 * np.pi / 180.0) * R * np.cos(lat1 * np.pi / 180.0)).value*1000, 3)
x_E34 = np.round(((delta_lon24 * np.pi / 180.0) * R * np.cos(lat1 * np.pi / 180.0)).value*1000, 3)
x_E43 = np.round(((delta_lon42 * np.pi / 180.0) * R * np.cos(lat1 * np.pi / 180.0)).value*1000, 3)

x_N12 = np.round((delta_lat12 * R * np.pi / 180.0).value*1000, 3)
x_N21 = np.round((delta_lat21 * R * np.pi / 180.0).value*1000, 3)
x_N13 = np.round((delta_lat13 * R * np.pi / 180.0).value*1000, 3)
x_N31 = np.round((delta_lat31 * R * np.pi / 180.0).value*1000, 3)
x_N14 = np.round((delta_lat14 * R * np.pi / 180.0).value*1000, 3)
x_N41 = np.round((delta_lat41 * R * np.pi / 180.0).value*1000, 3)
x_N23 = np.round((delta_lat23 * R * np.pi / 180.0).value*1000, 3)
x_N32 = np.round((delta_lat32 * R * np.pi / 180.0).value*1000, 3)
x_N24 = np.round((delta_lat24 * R * np.pi / 180.0).value*1000, 3)
x_N42 = np.round((delta_lat42 * R * np.pi / 180.0).value*1000, 3)
x_N34 = np.round((delta_lat34 * R * np.pi / 180.0).value*1000, 3)
x_N43 = np.round((delta_lat43 * R * np.pi / 180.0).value*1000, 3)

# Difference in height (x_up)
x_up12 = delta_height12.value
x_up21 = delta_height21.value
x_up13 = delta_height13.value
x_up31 = delta_height31.value
x_up14 = delta_height14.value
x_up41 = delta_height41.value
x_up23 = delta_height23.value
x_up32 = delta_height32.value
x_up24 = delta_height24.value
x_up42 = delta_height42.value
x_up34 = delta_height34.value
x_up43 = delta_height43.value


star_of_interest = "β Canis Majoris,Mirzam,0.019,52.632,0.77,1.98,1.75,28000.0,6.378,-16.044,06h 22m 42.0s,-17° 57′ 21″,0.69,0.57,0.52,1.6368e-05,7.967e-06,8.06e-06"

star_of_interest= "ε Orionis,Alnilam,-0.002,-500.0,0.47,1.7,1.51,30000.0,5.604,-0.798,05h 36m 12.8s,-01° 12′ 07″,0.76,0.63,0.56,2.1577e-05,1.0311e-05,1.0053e-05"
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

U12 = []
V12 = []
W12 = []

U21 = []
V21 = []
W21 = []

U13 = []
V13 = []
W13 = []

U31 = []
V31 = []
W31 = []

U14 = []
V14 = []
W14 = []

U41 = []
V41 = []
W41 = []

U23 = []
V23 = []
W23 = []

U32 = []
V32 = []
W32 = []

U24 = []
V24 = []
W24 = []

U42 = []
V42 = []
W42 = []

U43 = []
V43 = []
W43 = []

U34 = []
V34 = []
W34 = []

for time in tqdm(times):
    HA_value = RA_2_HA(given_ra_decimal, time.jd)

    matrices12 = R_x(given_dec_decimal * u.deg).dot(R_y(HA_value * u.deg)).dot(R_x(-lat1 * u.deg))
    h_plane12 = np.array([[x_E12], [x_N12], [x_up12]])  # x_N12 -> x_E12, x_E12 -> x_N12
    UVW_plane12 = matrices12.dot(h_plane12)

    U12.append(UVW_plane12[0][0])
    V12.append(UVW_plane12[1][0])
    W12.append(UVW_plane12[2][0])

    matrices21 = R_x(given_dec_decimal * u.deg).dot(R_y(HA_value * u.deg)).dot(R_x(-lat2 * u.deg))
    h_plane21 = np.array([[x_E21], [x_N21], [x_up21]])  # x_N21 -> x_E21, x_E21 -> x_N21
    UVW_plane21 = matrices21.dot(h_plane21)

    U21.append(UVW_plane21[0][0])
    V21.append(UVW_plane21[1][0])
    W21.append(UVW_plane21[2][0])

    matrices13 = R_x(given_dec_decimal * u.deg).dot(R_y(HA_value * u.deg)).dot(R_x(-lat1 * u.deg))
    h_plane13 = np.array([[x_E13], [x_N13], [x_up13]])  # x_N13 -> x_E13, x_E13 -> x_N13
    UVW_plane13 = matrices13.dot(h_plane13)

    U13.append(UVW_plane13[0][0])
    V13.append(UVW_plane13[1][0])
    W13.append(UVW_plane13[2][0])

    matrices31 = R_x(given_dec_decimal * u.deg).dot(R_y(HA_value * u.deg)).dot(R_x(-lat3 * u.deg))
    h_plane31 = np.array([[x_E31], [x_N31], [x_up31]])  # x_N31 -> x_E31, x_E31 -> x_N31
    UVW_plane31 = matrices31.dot(h_plane31)

    U31.append(UVW_plane31[0][0])
    V31.append(UVW_plane31[1][0])
    W31.append(UVW_plane31[2][0])

    matrices14 = R_x(given_dec_decimal * u.deg).dot(R_y(HA_value * u.deg)).dot(R_x(-lat1 * u.deg))
    h_plane14 = np.array([[x_E14], [x_N14], [x_up14]])  # x_N14 -> x_E14, x_E14 -> x_N14
    UVW_plane14 = matrices14.dot(h_plane14)

    U14.append(UVW_plane14[0][0])
    V14.append(UVW_plane14[1][0])
    W14.append(UVW_plane14[2][0])

    matrices41 = R_x(given_dec_decimal * u.deg).dot(R_y(HA_value * u.deg)).dot(R_x(-lat4 * u.deg))
    h_plane41 = np.array([[x_E41], [x_N41], [x_up41]])  # x_N41 -> x_E41, x_E41 -> x_N41
    UVW_plane41 = matrices41.dot(h_plane41)

    U41.append(UVW_plane41[0][0])
    V41.append(UVW_plane41[1][0])
    W41.append(UVW_plane41[2][0])

    matrices23 = R_x(given_dec_decimal * u.deg).dot(R_y(HA_value * u.deg)).dot(R_x(-lat2 * u.deg))
    h_plane23 = np.array([[x_E23], [x_N23], [x_up23]])  # x_N23 -> x_E23, x_E23 -> x_N23
    UVW_plane23 = matrices23.dot(h_plane23)

    U23.append(UVW_plane23[0][0])
    V23.append(UVW_plane23[1][0])
    W23.append(UVW_plane23[2][0])

    matrices32 = R_x(given_dec_decimal * u.deg).dot(R_y(HA_value * u.deg)).dot(R_x(-lat3 * u.deg))
    h_plane32 = np.array([[x_E32], [x_N32], [x_up32]])  # x_N32 -> x_E32, x_E32 -> x_N32
    UVW_plane32 = matrices32.dot(h_plane32)

    U32.append(UVW_plane32[0][0])
    V32.append(UVW_plane32[1][0])
    W32.append(UVW_plane32[2][0])

    matrices34 = R_x(given_dec_decimal * u.deg).dot(R_y(HA_value * u.deg)).dot(R_x(-lat3 * u.deg))
    h_plane34 = np.array([[x_E34], [x_N34], [x_up34]])
    UVW_plane34 = matrices34.dot(h_plane34)

    U34.append(UVW_plane34[0][0])
    V34.append(UVW_plane34[1][0])
    W34.append(UVW_plane34[2][0])

    matrices43 = R_x(given_dec_decimal * u.deg).dot(R_y(HA_value * u.deg)).dot(R_x(-lat4 * u.deg))
    h_plane43 = np.array([[x_E43], [x_N43], [x_up43]])
    UVW_plane43 = matrices43.dot(h_plane43)

    U43.append(UVW_plane43[0][0])
    V43.append(UVW_plane43[1][0])
    W43.append(UVW_plane43[2][0])

    matrices24 = R_x(given_dec_decimal * u.deg).dot(R_y(HA_value * u.deg)).dot(R_x(-lat2 * u.deg))
    h_plane24 = np.array([[x_E24], [x_N24], [x_up24]])
    UVW_plane24 = matrices24.dot(h_plane24)

    U24.append(UVW_plane24[0][0])
    V24.append(UVW_plane24[1][0])
    W24.append(UVW_plane24[2][0])

    matrices42 = R_x(given_dec_decimal * u.deg).dot(R_y(HA_value * u.deg)).dot(R_x(-lat4 * u.deg))
    h_plane42 = np.array([[x_E42], [x_N42], [x_up42]])
    UVW_plane42 = matrices42.dot(h_plane42)

    U42.append(UVW_plane42[0][0])
    V42.append(UVW_plane42[1][0])
    W42.append(UVW_plane42[2][0])

wavelength = 416e-9
Mlambda = wavelength*1000000
print(Mlambda)
plt.plot(np.array(U14) / Mlambda, np.array(V14) / Mlambda, '.', label="T1T4", color='y')
plt.plot(np.array(U41) / Mlambda, np.array(V41) / Mlambda, '.', color='y')
plt.plot(np.array(U12) / Mlambda, np.array(V12) / Mlambda, '.', label="T1T2", color='r')
plt.plot(np.array(U21) / Mlambda, np.array(V21) / Mlambda, '.', color='r')
plt.plot(np.array(U13) / Mlambda, np.array(V13) / Mlambda, '.', label="T1T3", color='b')
plt.plot(np.array(U31) / Mlambda, np.array(V31) / Mlambda, '.', color='b')
plt.plot(np.array(U23) / Mlambda, np.array(V23) / Mlambda, '.', label="T2T3", color='g')
plt.plot(np.array(U32) / Mlambda, np.array(V32) / Mlambda, '.', color='g')
plt.plot(np.array(U24) / Mlambda, np.array(V24) / Mlambda, '.', label="T2T4", color='m')
plt.plot(np.array(U42) / Mlambda, np.array(V42) / Mlambda, '.', color='m')
plt.plot(np.array(U34) / Mlambda, np.array(V34) / Mlambda, '.', label="T3T4", color='k')
plt.plot(np.array(U43) / Mlambda, np.array(V43) / Mlambda, '.', color='k')

plt.xlim(-450, 450)
plt.ylim(-350, 350)
plt.xlabel("u ($M\lambda$)")
plt.ylabel("v ($M\lambda$)")
plt.legend()
plt.title(BayerF+" "+date_str)
plt.show()