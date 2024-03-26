"""
Created on: 28.02.2023
Created by: Lucijana Stanic

"""
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

from brightstar_input import lat_deg1, lon_deg1, date_str, n_brightest_stars
from brightstar_functions import dms_to_decimal

import csv

from itertools import zip_longest

# - * - coding: utf - 8 - * -

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------


# Initialize a dictionary to store data from each column
data = {}

# Convert coordinates from deg to dec
lat_dec1 = dms_to_decimal(lat_deg1)
lon_dec1 = dms_to_decimal(lon_deg1)


# Open the CSV file
with open('stars_visible_' + date_str + '.csv', newline='') as csvfile:
    reader = csv.reader(csvfile)
    header = next(reader)  # Read the header row

    # Create an empty list for each column
    for column_name in header:
        data[column_name] = []

    # Iterate over each row in the CSV file
    for row in reader:
        # Iterate over each item in the row and append it to the respective column list
        for idx, item in enumerate(row):
            column_name = header[idx]  # Get the corresponding column name
            data[column_name].append(item)

number_of_stars = n_brightest_stars

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------


magnitude = input("Do you want to specify the apparent magnitude (V) of the star? (y/n) ")
if magnitude == "y":
    magnitude_threshold = float(input("Apparent magnitude threshold <= "))

    # Initialize extracted_data1 dictionary
    extracted_data1 = {}
    for key in data.keys():
        extracted_data1[key] = []

    # Loop over each star's Vmag and check against the threshold
    for star, vmag in enumerate(data["Vmag"]):
        if float(vmag) <= magnitude_threshold:
            # Append corresponding entries for the star to extracted_data1
            for key in data.keys():
                extracted_data1[key].append(data[key][star])

    print("Out of the ", number_of_stars, " analysed, ", len(extracted_data1[next(iter(data.keys()))]), " have "
                                                                                                        "an apparent "
                                                                                                        "magnitude <= "
                                                                                                        "",
          magnitude_threshold)
    list_2_print = input("Print list of those? (y/n)")
    if list_2_print == "y":
        print(extracted_data1)

    list_2_save = input("Save as .csv? (y/n)")
    if list_2_save == "y":
        output_file = 'stars_visible_bright_' + date_str + '.csv'

        # Write the extracted data to the CSV file
        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)

            # Write the header row using the keys of the extracted_data1 dictionary
            writer.writerow(extracted_data1.keys())

            # Write the data rows
            # Determine the number of rows based on the length of one of the lists in extracted_data1
            num_rows = len(next(iter(extracted_data1.values())))

            for i in range(num_rows):
                # Get the data for each column and write it to the CSV file
                row_data = [extracted_data1[key][i] for key in extracted_data1.keys()]
                writer.writerow(row_data)
        print("Saved as", output_file)
    extracted_data = extracted_data1


# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

# In case one would like to extract stars by estimated diameter
diameter_choice = input("Do you want to specify the diameter in mas of the star? (y/n) ")

if diameter_choice == "y":
    diameter_threshold = float(input("Diameter in mas >= "))

    # Initialize extracted_data2 dictionary
    extracted_data2 = {}
    for key in data.keys():
        extracted_data2[key] = []

    # Loop over each star's Vmag and check against the threshold
    for idx, (diam_V, diam_B, diam_U) in enumerate(
            zip_longest(extracted_data["Diameter_V"], extracted_data["Diameter_B"],
                        extracted_data["Diameter_U"])):
        # Convert diameters to float, handling empty strings
        diam_V = float(diam_V) if diam_V else None
        diam_B = float(diam_B) if diam_B else None
        diam_U = float(diam_U) if diam_U else None

        range_num = 0.1

        # Check if any of the diameters are within the threshold range_num
        if (diam_V is not None and abs(diam_V - diameter_threshold) <= range_num) or \
                (diam_B is not None and abs(diam_B - diameter_threshold) <= range_num) or \
                (diam_U is not None and abs(diam_U - diameter_threshold) <= range_num):
            # Append corresponding entries for the star to extracted_data2
            for key in data.keys():
                extracted_data2[key].append(extracted_data[key][idx])

    if magnitude == "y":
        print("Out of the ", number_of_stars, " analysed, ", len(extracted_data2[next(iter(data.keys()))]),
              "  have an apparent magnitude <= ", magnitude_threshold,
              "and a diameter of around ", diameter_threshold, " mas.")
    else:
        print("Out of the ", number_of_stars, " analysed, ", len(extracted_data2[next(iter(data.keys()))]),
              "  have a diameter of around ", diameter_threshold, " mas.")
    list_3_print = input("Print list of those? (y/n)")
    if list_3_print == "y":
        print(extracted_data2)

    list_3_save = input("Save as .csv? (y/n)")
    if list_3_save == "y":
        if magnitude == "y":
            output_file = 'stars_visible_bright_big_' + date_str + '.csv'
        else:
            output_file = 'stars_visible_big_' + date_str + '.csv'

        # Write the extracted data to the CSV file
        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)

            # Write the header row using the keys of the extracted_data2 dictionary
            writer.writerow(extracted_data2.keys())

            # Write the data rows
            # Determine the number of rows based on the length of one of the lists in extracted_data2
            num_rows = len(next(iter(extracted_data2.values())))

            for i in range(num_rows):
                # Get the data for each column and write it to the CSV file
                row_data = [extracted_data2[key][i] for key in extracted_data2.keys()]
                writer.writerow(row_data)
        print("Saved as", output_file)
    extracted_data = extracted_data1

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

# In case one would like to extract stars by certain Spectral photon flux Phi?
Phi_choice = input("Do you want to specify the spectral photon flux in photons m^-2 s^-1 Hz^-1 of the star? (y/n) ")

if Phi_choice == "y":
    Phi_threshold = float(input("Spectral photon flux Phi in photons m^-2 s^-1 Hz^-1 >= "))

    # Initialize extracted_data3 dictionary
    extracted_data3 = {}
    for key in data.keys():
        extracted_data3[key] = []

    # Loop over each star's Vmag and check against the threshold
    for idx, (Phi_V, Phi_B, Phi_U) in enumerate(
            zip_longest(extracted_data["Phi_V"], extracted_data["Phi_B"],
                        extracted_data["Phi_U"])):
        # Convert Phis to float, handling empty strings
        Phi_V = float(Phi_V) if Phi_V else None
        Phi_B = float(Phi_B) if Phi_B else None
        Phi_U = float(Phi_U) if Phi_U else None



        # Check if any of the Phis are within the threshold range_num
        if (Phi_V is not None and Phi_V >= Phi_threshold) or \
                (Phi_B is not None and Phi_B >= Phi_threshold) or \
                (Phi_U is not None and Phi_U >= Phi_threshold):
            # Append corresponding entries for the star to extracted_data3
            for key in data.keys():
                extracted_data3[key].append(extracted_data[key][idx])

    if magnitude == "y" and diameter_choice == "y":
        print("Out of the ", number_of_stars, " analysed, ", len(extracted_data3[next(iter(data.keys()))]),
              "  have an apparent magnitude <= ", magnitude_threshold,
              ", a diameter of around ", diameter_threshold, " mas and a Phi of around ", Phi_threshold ," photons m^-2 s^-1 Hz^-1.")

    if magnitude == "n" and diameter_choice == "y":
        print("Out of the ", number_of_stars, " analysed, ", len(extracted_data3[next(iter(data.keys()))]),
              "  have a diameter of around ", diameter_threshold, " mas and a Phi of around ", diameter_threshold ,"photons m^-2 s^-1 Hz^-1.")

    if magnitude == "y" and diameter_choice == "n":
        print("Out of the ", number_of_stars, " analysed, ", len(extracted_data3[next(iter(data.keys()))]),
              "  have an apparent magnitude <= ", magnitude_threshold,
              " and a Phi of around ", Phi_threshold ,"photons m^-2 s^-1 Hz^-1.")
    else:
        print("Out of the ", number_of_stars, " analysed, ", len(extracted_data3[next(iter(data.keys()))]),
              "  have a Phi of around ", Phi_threshold, " photons m^-2 s^-1 Hz^-1.")

    list_4_print = input("Print list of those? (y/n)")
    if list_4_print == "y":
        print(extracted_data3)

    list_4_save = input("Save as .csv? (y/n)")
    if list_4_save == "y":
        if magnitude == "y" and diameter_choice == "y":
            output_file = 'stars_visible_bright_big_highPhi' + date_str + '.csv'

        if magnitude == "n" and diameter_choice == "y":
            output_file = 'stars_visible_big_highPhi' + date_str + '.csv'

        if magnitude == "y" and diameter_choice == "n":
            output_file = 'stars_visible_bright_highPhi' + date_str + '.csv'

        else:
            output_file = 'stars_visible_highPhi_' + date_str + '.csv'

        # Write the extracted data to the CSV file
        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)

            # Write the header row using the keys of the extracted_data dictionary
            writer.writerow(extracted_data3.keys())

            # Write the data rows
            # Determine the number of rows based on the length of one of the lists in extracted_data
            num_rows = len(next(iter(extracted_data3.values())))

            for i in range(num_rows):
                # Get the data for each column and write it to the CSV file
                row_data = [extracted_data3[key][i] for key in extracted_data3.keys()]
                writer.writerow(row_data)
        print("Saved as", output_file)
    extracted_data = extracted_data3
last_output_file = output_file