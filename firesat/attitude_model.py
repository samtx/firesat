import numpy as np
import firesat.atmosphere as atmos

def attitude(x=None, y=None, var_info=None, fidelity=0, **kwargs):
    """Attitude control model to compute the torques necessary to counteract
    moments on the satellite based on perturbations. Then computes the power
    required to the necessary power required to apply torque to reaction
    control wheels.

    Parameters
    ----------
    x : np.ndarray (6, n)
        Input design vars
        x[0] = H
        x[1] = F_s
        x[2] = L_sp
        x[3] = q
        x[4] = L_a
        x[5] = C_d

    y : np.ndarray (3 or 5, n)
        Input coupling vars
        y[0] = v
        y[1] = dt_orbit
        y[2] = theta_slew
        y[3] = I_max  (optional)
        y[4] = I_min  (optional)

    var_info : dict
        Dictionary containing fixed parameters for problem

    Returns
    -------
    q : np.ndarray (2, n)
        Quantities of interest and output coupling variables
        q[0] = tau_tot
        q[1] = PACS

    Notes
    -----
    Idea for hifidelity version: calculate the atmospheric density, rho, based on elevation H
    """

    # Setup fixed input variables
    if not var_info:
        var_info = setup()

    RE = var_info["RE"]
    mu = var_info["mu"]
    c = var_info["c"]
    A_s = var_info["A_s"]
    sun_i = var_info["sun_i"]
    dt_slew = var_info["dt_slew"]
    M = var_info["M"]
    rho = var_info["rho"]
    A = var_info["A"]
    n = var_info["n"]
    omega_max = var_info["omega_max"]
    P_hold = var_info["P_hold"]
    theta = var_info["theta"]
    RD = var_info["RD"]

    # Unpack design variables
    if (x is not None) and (x.shape[0] == 6):
        H, F_s, L_sp, q, L_a, C_d = x
    else:
        H = kwargs.get("H", var_info["H_mean"])
        F_s = kwargs.get("F_s", var_info["F_s_mean"])
        L_sp = kwargs.get("L_sp", var_info["L_sp_mean"])
        q = kwargs.get("q", var_info["q_mean"])
        L_a = kwargs.get("L_a", var_info["L_a_mean"])
        C_d = kwargs.get("C_d", var_info["C_d_mean"])

    # Unpack coupling variables
    if (y is not None) and (y.shape[0] >= 3):
        v, dt_orbit, theta_slew = y[[0, 1, 2]]
        if y.shape[0] > 3:
            I_max, I_min = y[[3, 4]]
        else:
            I_max = var_info["I_max"]
            I_min = var_info["I_min"]
    else:
        v = kwargs.get("v", 0.5 * (var_info["v_min"] + var_info["v_max"]))
        dt_orbit = kwargs.get(
            "dt_orbit", 0.5 * (var_info["dt_orbit_min"] + var_info["dt_orbit_max"])
        )
        theta_slew = kwargs.get(
            "theta_slew",
            0.5 * (var_info["theta_slew_min"] + var_info["theta_slew_max"]),
        )
        I_max = kwargs.get("I_max", var_info["I_max"])
        I_min = kwargs.get("I_min", var_info["I_min"])

    # Compute output quantitites

    if fidelity == 0:
        # Low fidelity computation
        tau_slew = 4 * theta_slew / (dt_slew) ** 2 * I_max
        tau_g = 3 * mu / (2 * ((RE + H) ** 3)) * abs(I_max - I_min) * np.sin(2 * theta)
        tau_sp = L_sp * F_s / c * A_s * (1 + q) * np.cos(sun_i)
        tau_m = 2 * M * RD / ((RE + H) ** 3)
        tau_a = 0.5 * L_a * rho * C_d * A * v ** 2
        tau_dist = np.sqrt(tau_g ** 2 + tau_sp ** 2 + tau_m ** 2 + tau_a ** 2)
        tau_tot = np.maximum(tau_slew, tau_dist)
        # tau_tot    = tau_slew + tau_dist
        PACS = tau_tot * omega_max + n * P_hold

    if fidelity > 0:
        # Medium fidelity computation
        tau_slew = 4 * theta_slew / (dt_slew) ** 2 * I_max
        tau_g = 3 * mu / (2 * ((RE + H) ** 3)) * abs(I_max - I_min) * np.sin(2 * theta)
        tau_sp = L_sp * F_s / c * A_s * (1 + q) * np.cos(sun_i)
        tau_m = 2 * M * RD / ((RE + H) ** 3)

        # use exponential atmospheric density model
        rho = atmos.exponential_density_model(H/1000)
        tau_a = 0.5 * L_a * rho * C_d * A * v ** 2


        tau_dist = np.sqrt(tau_g ** 2 + tau_sp ** 2 + tau_m ** 2 + tau_a ** 2)
        tau_tot = np.maximum(tau_slew, tau_dist)
        # tau_tot    = tau_slew + tau_dist
        PACS = tau_tot * omega_max + n * P_hold


    # Assemble Attitude Control outputs
    n = x.shape[1]  # number of samples
    q = np.zeros((2, n))
    q[0] = tau_tot
    q[1] = PACS
    return q
