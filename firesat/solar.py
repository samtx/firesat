# Functions to determine solar-satellite vector

import numpy as np
from numpy.linalg import norm
from firesat.constants import DEG2RAD, AU_KM, R_EARTH
import math


def sun_sat_angle(rsat, rsun):
    """Compute the sun-satellite angle
    Args:
        rsat : float (n, 3), satellite position vector in ECI [km]
        rsun : float (n, 3), sun position vector in ECI [km]
    Outputs:
        zeta : float (n), angle in radians
    References:
        Vallado, p. 912, Alg. 74
    """
    # Sun-satellite angle
    rsat = np.atleast_2d(rsat)
    rsun = np.atleast_2d(rsun)
    n = rsat.shape[0]
    sinzeta = np.empty(n)
    # for i in range(n):
    #     crosspdt = np.cross(rsun[i], rsat[i])
    #     sinzeta[i] = norm(crosspdt)/(norm(rsun[i])*norm(rsat[i]))
    crosspdt = np.cross(rsun, rsat, axis=1)
    sinzeta = norm(crosspdt, axis=1)/(norm(rsun, axis=1)*norm(rsat,axis=1))
    return np.arcsin(sinzeta)


def sun_sat_orthogonal_distance(rsat, zeta):
    """
    Args:
        rsat : float (n, 3)
        zeta : float (n)
    Output:
        float (n) : distance in km
    """
    tmp1 = norm(rsat, axis=1)
    tmp2 = np.cos(zeta-math.pi*0.5)
    return tmp1 * tmp2


def is_sat_illuminated(rsat, rsun):
    """Determine if satellite is illuminated by sun
    Args:
        rsat : float (n, 3)
        rsun : float (n, 3)
    Output:
        vis : bool (n)
    """
    zeta = sun_sat_angle(rsat, rsun)
    dist = sun_sat_orthogonal_distance(rsat, zeta)
    return dist > R_EARTH


def sun_pos(jdt):
    """Compute the Sun position vector from julian date
    Args:
        jdt : float (n), array or scalar of julian dates
    Output:
        r : float (n, 3), position vector of sun in ECI coordinates [km]
    References:
        Vallado, p. 279, Alg. 29
        Vallado software, AST2BODY.FOR, subroutine SUN
    """
    jdt = np.atleast_1d(jdt)
    t_ut1 = (jdt - 2451545.0)/36525
    t_tdb = t_ut1
    lmda_Msun = (280.4606184 + 36000.77005361*t_tdb) % 360
    # M_sun = (357.5291092 + 35999.05034*t_tdb) % 360
    M_sun = (357.5277233 + 35999.05034*t_tdb) % 360
    lmda_eclp = lmda_Msun + 1.914666471*np.sin(M_sun*DEG2RAD)
    lmda_eclp += 0.019994643*np.sin(2*M_sun*DEG2RAD)
    r_sun_mag = 1.000140612 - 0.016708617*np.cos(M_sun*DEG2RAD)
    r_sun_mag -= 0.000139589*np.cos(2*M_sun*DEG2RAD)
    eps = 23.439291 - 0.0130042*t_tdb
    coslmda = np.cos(lmda_eclp*DEG2RAD)
    sinlmda = np.sin(lmda_eclp*DEG2RAD)
    coseps = np.cos(eps*DEG2RAD)
    sineps = np.sin(eps*DEG2RAD)
    if jdt.size > 1:
        r = np.empty((3,jdt.size))
    else:
        r = np.empty(3)
    r[0] = r_sun_mag * coslmda
    r[1] = r_sun_mag * coseps * sinlmda
    r[2] = r_sun_mag * sineps * sinlmda
    r *= AU_KM
    return r.T
