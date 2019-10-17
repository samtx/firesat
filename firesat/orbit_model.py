import numpy as np
from firesat import sgp4, solar, timefn
import firesat.constants as cst

def orbit(x, var_info, fidelity=0):
    """Calculate the orbit of the satellite based on the height.

    Parameters
    ----------
    x : np.ndarray (2, n)
        Input design vars
        x[0] = H
        x[1] = phi (optional)

    var_info : dict
        Dictionary containing fixed parameters for problem

    Returns
    -------
    q : np.ndarray (4, n)
        Quantities of interest and output coupling variables of orbit model
        q[0] = v
        q[1] = dt_orbit
        q[2] = dt_eclipse
        q[3] = theta_slew

    Notes
    -----
    Idea for hifidelity version: calculate the eclipse time based on
    integrating an circular/elliptical orbit the given elevation. Can change
    orbital parameters to have additional input variables. Can calculate if the
    satellite is in eclipse based on an ellipsoidal Earth model.

    Based on real satellite Terra, TLE below
    TERRA
    1 25994U 99068A   19289.08873013  .00000005  00000-0  11240-4 0  9997
    2 25994  98.1961   1.7962 0001473  92.0911 268.0459 14.57115279 54585
    """
    # Setup fixed variables
    mu = var_info["mu"]
    RE = var_info["RE"]

    if x.ndim == 1:
        H = x
        phi = var_info["phi"]
        n = x.size  # number of samples
    else:
        H, phi = x
        n = x.shape[1]

    # Compute Orbit Subsystem Outputs
    if fidelity == 0:
        # Low fidelity model
        v = np.sqrt(mu / (RE + H))
        dt_orbit = 2 * np.pi * (RE + H) / v
        dt_eclipse = dt_orbit / np.pi * np.arcsin(RE / (RE + H))
        theta_slew = np.arctan(np.sin(phi / RE) / (1 - np.cos(phi / RE) + H / RE))
    else:
        # Compute high fidelity model
        # mean motion [rev/min]
        no = np.sqrt(mu / ((RE + H) ** 3)) * 60
        # propagate using sgp4
        t = np.linspace(0, 1440, 1441)
        tsize = t.size
        satrec = sgp4.Satellite()
        v_all = np.empty(n)
        dt_orbit_all = np.empty(n)
        dt_eclipse_all = np.empty(n)
        theta_slew_all = np.empty(n)
        for i in range(n):
            satrec.no = no[i]
            sgp4.sgp4init(
                satrec.whichconst,
                False,  # afspc_mode = False
                satrec.satnum,
                satrec.epoch,
                satrec.bstar,
                satrec.ecco,
                satrec.argpo,
                satrec.inclo,
                satrec.mo,
                satrec.no,
                satrec.nodeo,
                satrec,
            )
            r = np.empty((tsize, 3))
            v = np.empty((tsize, 3))
            for j in range(t.size):
                ri, vi = sgp4.sgp4(satrec, t[j])
                r[j] = ri
                v[j] = vi

            # compute velocity
            vnorms = np.linalg.norm(v, axis=1)
            v_all[i] = np.mean(vnorms) * 1000  # km/s --> m/s

            # compute orbit time
            # find indecies where y-axis value goes from negative to positive, then linearly interpolate for times
            # this function finds the sign change, and where r2[i+1]>r2[i]
            idx = np.nonzero((r[:-1, 1] * r[1:, 1] < 0) & (r[:-1, 1] < r[1:, 1]))[0]
            t_r20 = np.diff(t)[idx] / np.diff(r[:, 1])[idx] * (-r[idx, 1]) + t[idx]
            dt_orbit = np.diff(t_r20) * 60  # min --> sec
            dt_orbit_avg = np.mean(dt_orbit)
            dt_orbit_all[i] = dt_orbit_avg

            # compute eclipse time
            epoch = satrec.epoch
            jdt = timefn.jdt_tsince(epoch+cst.J2000, t)
            rsun = solar.sun_pos(jdt)
            # count the number of visible time steps per orbit
            ecl = np.empty(len(idx)-1)
            for j in range(len(idx)-1):
                idx1, idx2 = idx[j], idx[j+1]
                vis = solar.is_sat_illuminated(r[idx1:idx2], rsun[idx1:idx2])
                ecl[j] = 1 - np.sum(vis.astype(int))/(idx2 - idx1)
            dt_eclipse_all[i] = np.mean(ecl) * dt_orbit_avg

            # compute slewing angle
            rnorm = np.linalg.norm(r, axis=1)
            H_i = rnorm*1000 - RE  # altitude at time t
            theta_slew = np.arctan(np.sin(phi[i] / RE) / (1 - np.cos(phi[i] / RE) + H_i / RE))
            theta_slew_all[i] = np.mean(theta_slew)
            # print(f'i = {i}')
        v = v_all
        dt_orbit = dt_orbit_all
        dt_eclipse = dt_eclipse_all
        theta_slew = theta_slew_all

    # Assemble Orbit outputs
    q = np.zeros((4, n))
    q[0] = v
    q[1] = dt_orbit
    q[2] = dt_eclipse
    q[3] = theta_slew
    return q
