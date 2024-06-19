"""
Created on: 11.06.2023
Created by: Lucijana Stanic

"""
# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy.constants import c, k, h, pi
import matplotlib.colors as mcolors

SMALL_SIZE = 8
MEDIUM_SIZE = 13
BIGGER_SIZE = 18

plt.rc('font', size=BIGGER_SIZE)  # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=BIGGER_SIZE)  # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

np.random.seed(3)

# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------


# Initialize a dictionary to store data from each column
data = {}
plt.style.use('dark_background')


# Open the CSV file
with open('10000stars_data.csv', newline='') as csvfile:
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


def convert_strings_to_floats(input_array):
    output_array = []
    for element in input_array:
        converted_float = float(element)
        output_array.append(converted_float)
    return output_array


Phi_V = np.array(data['Phi_V'])
Phi_V = convert_strings_to_floats(Phi_V)

diameter_V = np.array(data['Diameter_V'])
diameter_V = convert_strings_to_floats(diameter_V)

temps = np.array(data['Temp'])
temps = convert_strings_to_floats(temps)

inverse_diameter = []
for i in range(len(diameter_V)):
    inverse_diameter.append(1 / diameter_V[i])

# Plot using the colormap based on temperatures
sc = plt.scatter(inverse_diameter, Phi_V, c=temps, cmap=bb_cmap, marker='.')
plt.colorbar(sc, label='Temperature (K)')
plt.yscale('log')
plt.xlabel('1/θ')
plt.ylabel('Φ')
plt.title('Φ vs θ')
plt.ylim(1e-12, 1e-2)
plt.show()


def Phi_vs_theta(wavelength, T, theta):
    top = theta ** 2 * pi ** 3
    bottom = 36 * 60 ** 6 * wavelength ** 2 * np.exp(h * c / (wavelength * k * T))
    return top / bottom


temperatures = [1000, 2000, 3000, 5000, 10000]
wavelengths = [364 * 1e-9, 442 * 1e-9, 540 * 1e-9]
wavelengths = [540 * 1e-9]

thetas = np.linspace(0.00015, 0.05, 100)

# Define colors and linestyles for temperatures and wavelengths
colors_rgb = ['#ff3800', '#ff8912', '#ffb46b', '#ffe4ce', '#ccdbff']
linestyles = ['-', '--', ':']


for i, T in enumerate(temperatures):
    for j, wavelength in enumerate(wavelengths):
        color = colors_rgb[i]
        linestyle = linestyles[j]
        label = f'λ={wavelength * 1e9:.0f}nm, T={T}K'
        plt.semilogy(1 / thetas, Phi_vs_theta(wavelength, T, thetas), color=color, linestyle=linestyle, label=label)

# Add labels and legend
plt.xlabel('θ (arcseconds)')
plt.ylabel('Φ')
plt.title('Φ vs θ')
#plt.legend()
plt.grid(True)
plt.ylim(1e-12, 1e-2)
plt.show()
