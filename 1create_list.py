"""
Created on: 28.02.2023
Created by: Lucijana Stanic

This code retrieves the most important properties from the "bsc5-all.json" file for our purposes and saves them in
"stars_data.csv"
"""
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

import json
import csv
from brightstar_functions import process_star
from brightstar_input import n_brightest_stars

# - * - coding: utf - 8 - * -

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------


lambda_U = 364 * 10 ** (-9)
lambda_V = 540 * 10 ** (-9)
lambda_B = 442 * 10 ** (-9)

# Load the JSON data
with open('bsc5-all.json', 'r') as file:
    stars_data = json.load(file)

# Sort stars_data based on Vmag (magnitude) in descending order
brightest_stars = sorted(stars_data, key=lambda x: x['Vmag'])[:n_brightest_stars]

# Process stars data
data = [process_star(star) for star in brightest_stars]

# Write data to CSV file
csv_columns = ["BayerF", "Common", "Parallax", "Distance", "Umag", "Vmag", "Bmag", "Temp", "RA_decimal", "Dec_decimal",
               "RA", "Dec", "Diameter_U", "Diameter_V", "Diameter_B"]
csv_file = str(n_brightest_stars)+"stars_data.csv"
with open(csv_file, 'w') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
    writer.writeheader()
    writer.writerows(data)

print("CSV file with {} brightest stars generated successfully.".format(n_brightest_stars))