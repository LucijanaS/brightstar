import numpy as np
from scipy.optimize import brentq
from scipy.special import j1


class Orbit():

    def __init__(self, e, I, Omega, omega, tph):
        self.e, self.tph = e, tph,
        pi, cos, sin = np.pi, np.cos, np.sin
        x = I * pi / 180
        self.cosI = cos(x)
        x = Omega
        self.cosOm, self.sinOm = cos(x), sin(x)
        x = omega
        self.cosom, self.sinom = cos(x), sin(x)

    def pos(self, tph):
        ''' Keplerian x,y,z '''
        e = self.e
        # First solve Kepler's equation
        # tph is orbital phase (aka mean anomaly)
        psi = brentq(lambda psi: psi - e * np.sin(psi) - tph, 0, 2 * np.pi)
        # then find x,y in orbital plane
        ec = (1 - e * e) ** .5
        x, y = np.cos(psi) - e, ec * np.sin(psi)
        # next rotate by omega
        cs, sn = self.cosom, self.sinom
        x, y = x * cs - y * sn, x * sn + y * cs
        # then inclination
        x, y = x, y * self.cosI
        # and finally rotate by Omega
        cs, sn = self.cosOm, self.sinOm
        x, y = x * cs - y * sn, x * sn + y * cs
        return x, y

    def binarypos(self, tph):
        q = 1
        a1 = q / (1 + q)
        a2 = -1 / (1 + q)
        tph = tph % (2 * np.pi)
        xs, ys = self.pos(tph)
        return [a1 * xs, a2 * xs], [a1 * ys, a2 * ys]


def visibility_binaries(u, v, theta_1, theta_2, x_1, y_1, x_2, y_2, lambda_):
    """
    Computes the squared visibility |V_12|^2 for a binary star system
    as observed using interferometry.

    Parameters:
    -----------
    u, v : 2D arrays (spatial frequency grids)
        Spatial frequency components corresponding to the projected baseline (units of inverse meters, [1/m]).
    theta_1, theta_2 : float
        Angular diameters (in radians) of the two binary stars.
    x_1, y_1 : float
        Position (in radians) of the first binary star, relative to the system's center of mass.
    x_2, y_2 : float
        Position (in radians) of the second binary star, relative to the system's center of mass.
    lambda_ : float
        Wavelength of observation (in meters).

    lambda_ : float
        Wavelength of observation (in meters).

    Returns:
    --------
    V_2 : float or array-like
        Squared visibility (|V_12|^2), which represents the coherence of the light from the binary stars.
    """

    input_1 = np.pi * np.sqrt(u ** 2 + v ** 2) * theta_1 / lambda_
    V_ud_1 = (2 * j1(input_1) / input_1)

    input_2 = np.pi * np.sqrt(u ** 2 + v ** 2) * theta_2 / lambda_
    V_ud_2 = (2 * j1(input_2) / input_2)

    exp_1 = np.exp(((2 * np.pi * 1j) / lambda_) * (u * x_1 + v * y_1))
    exp_2 = np.exp(((2 * np.pi * 1j) / lambda_) * (u * x_2 + v * y_2))

    V_binaries = V_ud_1 * exp_1 + V_ud_2 * exp_2

    V_2 = abs(V_binaries) ** 2

    return V_2


def seeing_double(theta_A, theta_B, theta_a, e, I, omega, Omega, U, V, tph):
    """
    The visibility measured for a binary system given its sizes, semi-major axis, orbital parameters and U,V.
    All parameters are unitless.

    Parameters to be fitted:
    - theta_A (float): angular diameter of star A
    - theta_B (float): angular diameter of star B
    - theta_a (float): semi-major axis
    - e (float): Orbital eccentricity (0 <= e < 1).
    - I (float): Orbital inclination (0 <= I < 2pi).
    - omega (float): Argument of periapsis (0 <= I < 2pi).
    - Omega (float): Longitude of ascending node (0 <= I < 2pi).
    - phase (float): Phase of orbit (0 <= I < 2pi).

    Parameters known or input:
    - U (float): U coordinate in the UV plane (should be unitless, i.e. divide by wavelength).
    - V (float): V coordinate in the UV plane (should be unitless, i.e. divide by wavelength).

    Returns:
    - Visibility measured given the parameters.
    """

    binary_selected = Orbit(e, I, Omega, omega, tph)

    xp, yp = binary_selected.binarypos(tph)

    return visibility_binaries(U, V, theta_A, theta_B, xp[0] * theta_a, yp[0] * theta_a, xp[1] * theta_a,
                               yp[1] * theta_a, 1)
