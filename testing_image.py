import numpy as np

import matplotlib.pyplot as plt

from scipy.special import j1, j0


def intensity(b, theta, lambda_):
    input = np.pi * b * theta / lambda_
    B_1 = (2 * j1(input) / input) ** 2
    return B_1


def meshgrid(s, resolution):
    u = np.linspace(-s, s, resolution)
    a, b = np.meshgrid(u, u)
    m = (a * a + b * b) ** (1 / 2)
    return m


resolution = 300
diameter = 0.0001  # in arcseconds
diameter_in_rad = diameter * np.pi / (3600 * 180)  # in radians
diameter_mas = diameter * 1000  # in milliarcseconds
wavelength = 5.4e-7  # wavelength in meters
wavelength_nm = wavelength * 10 ** 9

# Create a grid of points
size_to_plot = 150

x = np.linspace(-size_to_plot, size_to_plot, resolution)
y = np.linspace(-size_to_plot, size_to_plot, resolution)
X, Y = np.meshgrid(x, y)
R = np.sqrt(X ** 2 + Y ** 2)

# Compute intensity for each point in the grid
intensity_values = intensity(R, diameter_in_rad, wavelength)

# Plot the intensity as a 2D image array radially
plt.figure(figsize=(8, 8))
plt.imshow(intensity_values, extent=(-size_to_plot, size_to_plot, -size_to_plot, size_to_plot), origin='lower')
plt.colorbar(label='Intensity')
plt.title('Resolution: {}, diameter: {}mas, wavelength: {}nm'.format(resolution, diameter_mas, wavelength_nm))
plt.show()
