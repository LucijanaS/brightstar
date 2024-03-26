"""
Created on: 20.03.2023
Created by: Lucijana Stanic

All the functions needed in this repository are defined in this

"""
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------


import json
import csv
from datetime import date, datetime
import pytz
from tqdm import tqdm
from geopy.geocoders import Nominatim
from timezonefinder import TimezoneFinder
from itertools import zip_longest

import numpy as np
from scipy.constants import c, h, k, pi
from scipy.special import j0, j1
import scipy.signal as signal

import matplotlib.pyplot as plt

import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

from scipy.spatial import ConvexHull

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

lambda_U = 364 * 10 ** (-9)
lambda_V = 540 * 10 ** (-9)
lambda_B = 442 * 10 ** (-9)


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

def Phi(mag, wavelength):
    if mag is not None:
        nu = c / wavelength
        return np.round(10 ** (-22.44 - mag / 2.5) / (2 * nu * h), 9)
    else:
        return None

def calculate_diameter(mag, wavelength, temp):
    if temp is not None and mag is not None:
        nu = c / wavelength
        S = (nu ** 2 / c ** 2) / np.exp((h * nu) / (k * temp))
        area_steradian = Phi(mag, wavelength) / S

        radius_radians = np.sqrt(area_steradian / (pi))
        diameter_ = (6 / pi) * 60 ** 3 * radius_radians
        diameter = np.round(diameter_ * 10 ** 3, 2)
    else:
        diameter = None
    return diameter


def process_star(star):
    BayerF = star.get("BayerF")
    common = star.get("Common")

    parallax_value = star.get("Parallax")
    parallax = float(parallax_value) if parallax_value is not None else None
    distance = round(1 / parallax, 3) if parallax_value and abs(parallax) > 0 else None
    Vmag = round(float(star.get("Vmag")), 3)

    BV_value = star.get("B-V")
    BV = float(BV_value) if BV_value is not None else None
    Bmag = round(BV + Vmag, 3) if BV is not None else None

    UB_value = star.get("U-B")
    UB = float(UB_value) if UB_value is not None and Bmag is not None else None
    Umag = round(UB + Bmag, 3) if UB is not None else None

    temp = round(float(star.get("K")), 3) if star.get("K") is not None else None

    star_ra_decimal, star_dec_decimal = convert_ra_dec(star["RA"], star["Dec"])
    star_ra_decimal = round(star_ra_decimal, 3)
    star_dec_decimal = round(star_dec_decimal, 3)

    diameter_U = calculate_diameter(Umag, lambda_U, temp)
    diameter_V = calculate_diameter(Vmag, lambda_V, temp)
    diameter_B = calculate_diameter(Bmag, lambda_B, temp)

    Phi_V = Phi(Vmag, lambda_V)
    Phi_B = Phi(Bmag, lambda_B)
    Phi_U = Phi(Umag, lambda_U)

    return {
        "BayerF": BayerF,
        "Common": common,
        "Parallax": parallax,
        "Distance": distance,
        "Umag": Umag,
        "Vmag": Vmag,
        "Bmag": Bmag,
        "Temp": temp,
        "RA_decimal": star_ra_decimal,
        "Dec_decimal": star_dec_decimal,
        "RA": star.get("RA"),
        "Dec": star.get("Dec"),
        "Diameter_U": diameter_U,
        "Diameter_V": diameter_V,
        "Diameter_B": diameter_B,
        "Phi_U": Phi_U,
        "Phi_V": Phi_V,
        "Phi_B": Phi_B
    }


def intensity(b, theta, lambda_):
    input = np.pi * b * theta / lambda_
    B_1 = (2 * j1(input) / input) ** 2
    return B_1