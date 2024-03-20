"""
Created on: 28.02.2023
Created by: Lucijana Stanic

This code retrieves the most important properties from the "bsc5-all.json" file for our purposes and saves them in
"stars_data.csv"
"""
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

import numpy as np
from scipy.constants import c, h, k, pi
import json
import csv

# - * - coding: utf - 8 - * -

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------


lambda_U = 364 * 10 ** (-9)
lambda_V = 540 * 10 ** (-9)
lambda_B = 442 * 10 ** (-9)

data = []


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


# Load the JSON data
with open('bsc5-all.json', 'r') as file:
    stars_data = json.load(file)

# Sort stars_data based on Vmag (magnitude) in descending order
n=1000
brightest_stars = sorted(stars_data, key=lambda x: x['Vmag'])[:n]

# Loop through the top 1000 stars
for star in brightest_stars:
    # Extract properties if available
    BayerF = star.get("BayerF")
    common = star.get("Common")

    # Handling parallax
    parallax_value = star.get("Parallax")
    if parallax_value is not None:
        parallax = float(parallax_value)
        if abs(parallax) > 0:
            distance = 1 / parallax
        else:
            distance = None
            parallax = None
    else:
        distance = None
        parallax = None

    Vmag = float(star.get("Vmag"))

    BV_value = star.get("B-V")
    if BV_value is not None:
        BV = float(BV_value)
        Bmag = BV + Vmag
    else:
        Bmag = None
        BV = None

    UB_value = star.get("U-B")
    if UB_value is not None and Bmag is not None:
        UB = float(UB_value)
        Umag = UB + Bmag
    else:
        Umag = None
        UB = None

    temp_value = star.get("K")
    if temp_value is not None:
        temp = float(temp_value)
    else:
        temp = None


    star_ra_decimal, star_dec_decimal = convert_ra_dec(star["RA"], star["Dec"])

    RA_str = star.get("RA")
    Dec_str = star.get("Dec")

    if temp is not None:
        nu_V = c / lambda_V
        Phi_V = 10 ** (-22.44 - Vmag / 2.5) / (2 * nu_V * h)
        S_V = (nu_V ** 2 / c ** 2) / np.exp((h * nu_V) / (k * temp))
        area_steradian_V = Phi_V / S_V

        radius_radians_V = np.sqrt(area_steradian_V / (pi))
        diameter_V_ = (6 / pi) * 60 ** 3 * radius_radians_V
        diameter_V = np.round(diameter_V_ * 10 ** 3, 2)
    else:
        diameter_V_ = None

    if Bmag is not None and temp is not None:
        nu_B = c / lambda_B
        Phi_B = 10 ** (-22.44 - Bmag / 2.5) / (2 * nu_B * h)
        S_B = (nu_B ** 2 / c ** 2) / np.exp((h * nu_B) / (k * temp))
        area_steradian_B = Phi_B / S_B

        radius_radians_B = np.sqrt(area_steradian_B / (pi))
        diameter_B_ = (6 / pi) * 60 ** 3 * radius_radians_B
        diameter_B = np.round(diameter_B_ * 10 ** 3, 2)
    else:
        diameter_B = None

    if Umag is not None and temp is not None:
        nu_U = c / lambda_U
        Phi_U = 10 ** (-22.44 - Umag / 2.5) / (2 * nu_U * h)
        S_U = (nu_U ** 2 / c ** 2) / np.exp((h * nu_U) / (k * temp))
        area_steradian_U = Phi_U / S_U

        radius_radians_U = np.sqrt(area_steradian_U / (pi))
        diameter_U_ = (6 / pi) * 60 ** 3 * radius_radians_U
        diameter_U = np.round(diameter_U_ * 10 ** 3, 2)
    else:
        diameter_U = None

    data.append({
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
        "RA": RA_str,
        "Dec": Dec_str,
        "Diameter_U": diameter_U,
        "Diameter_V": diameter_V,
        "Diameter_B": diameter_B,
    })

# Write data to CSV file
csv_columns = ["BayerF", "Common", "Parallax", "Distance", "Umag", "Vmag", "Bmag", "Temp", "RA_decimal", "Dec_decimal",
               "RA", "Dec", "Diameter_U", "Diameter_V", "Diameter_B"]
csv_file = str(n)+"stars_data.csv"
with open(csv_file, 'w') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
    writer.writeheader()
    for star_data in data:
        writer.writerow(star_data)

print("CSV file generated successfully.")
