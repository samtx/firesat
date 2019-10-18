# Atmospheric Models

import numpy as np

WERTZ_1978_ATMOS = np.array([
#   h_0 [km], rho_0 [kg/m3],    H [km]
    [    150,     2.070e-09,   22.523],
    [    180,     5.464e-10,   29.740],
    [    200,     2.789e-10,   37.105],
    [    250,     7.248e-11,   45.546],
    [    300,     2.418e-11,   53.628],
    [    350,     9.518e-12,   53.298],
    [    400,     3.725e-12,   58.515],
    [    450,     1.585e-12,   60.828],
    [    500,     6.967e-13,   63.822],
    [    600,     1.454e-13,   71.835],
    [    700,     3.614e-14,   88.667],
    [    800,     1.170e-14,  124.640],
    [    900,     5.245e-15,  181.050],
    [   1000,     3.019e-15,  268.000],
])

def exponential_density_model(h):
    """Exponential density model. Assumes a spherically symmetrical distribution
    of particles in which the density varies exponentially. Based on US standard
    atmosphere model for h = [0, 25km], CIRA-72 for h = [25km, 500km], and
    CIRA-72 with Tinf = 1000K for h >= 500km

    Args:
        h : float (n), height above the ellipsoid [km]
    Output:
        rho : float (n), atmospheric density at altitude [kg/m^3]

    Notes:
        Scale height, H, is the fractional change in density with height. Wertz (1978)
        shows that scale height is equal to (k)*(temperature)/(molecular weight)*(gravity)
        where k is the Boltzmann constant

    References:
        Vallado, p. 565, Eq. 8-33
    """
    idx = np.searchsorted(WERTZ_1978_ATMOS[:,0], h, side='right') - 1
    h0 = WERTZ_1978_ATMOS[idx, 0]
    rho0 = WERTZ_1978_ATMOS[idx, 1]
    H = WERTZ_1978_ATMOS[idx, 2]
    rho = rho0 * np.exp(-(h - h0)/H)
    return rho
