from numpy import sqrt, sin, cos, arctan, arcsin, pi
import numpy as np
from firesat import orbit, attitude, power


def setup():
    """Function returns a dict containing the values for the fixed inputs
    (deterministic values) and random input parameters (mean and
    standard deviation).

    The satellite problem is adapted from the reference:
    Sankararaman, S., and Mahadevan, S. Likelihood-based approach to
    multidisciplinary analysis under uncertainty.
    Journal of Mechanical Design, 2012.
    """

    # Define variables names
    rand_inputs = ["H", "Po", "F_s", "L_sp", "q", "L_a", "C_d", "phi"]
    fixed_inputs = {
        "RE",
        "mu",
        "c",
        "A_s",
        "sun_i",
        "dt_slew",
        "M",
        "rho",
        "A",
        "n",
        "omega_max",
        "P_hold",
        "I_d",
        "nu",
        "LT",
        "eps_deg",
        "r_lw",
        "n_sa",
        "rho_sa",
        "t",
        "D",
        "I_bodyx",
        "I_bodyy",
        "I_bodyz",
        "theta",
        "RD",
    }
    output_vars = {
        "v",
        "dt_orbit",
        "dt_eclipse",
        "theta_slew",
        "tau_slew",
        "tau_g",
        "tau_sp",
        "tau_m",
        "tau_a",
        "tau_dist",
        "tau_tot",
        "PACS",
        "P_tot",
        "T_e",
        "T_d",
        "P_e",
        "P_d",
        "P_sa",
        "PBOL",
        "PEOL",
        "A_sa",
        "L",
        "W",
        "m_sa",
        "I_sax",
        "I_say",
        "I_saz",
        "I_max",
        "I_min",
    }

    # Declare fixed quantities
    RE = 6378140
    mu = 3.986 * (10 ** 14)
    phi = 235
    c = 2.9989 * (10 ** 8)
    A_s = 13.85
    sun_i = 0
    dt_slew = 760
    M = 7.96 * (10 ** 15)
    rho = 5.1480 * (10 ** (-11))
    A = 13.85
    n = 3
    omega_max = 6000
    P_hold = 20
    I_d = 0.77
    nu = 0.22
    LT = 15
    eps_deg = 0.0375
    r_lw = 3
    n_sa = 3
    rho_sa = 700
    t = 0.005
    D = 2
    I_bodyx = 6200
    I_bodyy = 6200
    I_bodyz = 4700
    theta = 15
    RD = 5

    # Mean value for coupling variables
    I_max = 6612.9
    I_min = 5116.0

    # Declare random quantity info
    H_mean = 18000000
    H_std = 1000000

    phi_mean = phi
    phi_std = 10

    Po_mean = 1000
    Po_std = 50

    F_s_mean = 1400
    F_s_std = 20

    # theta_mean = 15
    # theta_std  = 1

    L_sp_mean = 2
    L_sp_std = 0.4

    q_mean = 0.5
    q_std = 0.1

    # RD_mean    = 5
    # RD_std     = 1

    L_a_mean = 2
    L_a_std = 0.4

    C_d_mean = 1
    C_d_std = 0.2

    # Coupling variable uniform distribution min/max

    # if H_mean = 1800000 (18e5)
    # PACS_min        = 30.3
    # PACS_max        = 1085.
    # dt_orbit_min    = 3340.
    # dt_orbit_max    = 74229.
    # dt_eclipse_min  = 1044.
    # dt_eclipse_max  = 5215.
    # theta_slew_min  = 0.00000462
    # theta_slew_max  = 0.0002839
    # v_min           = 1849.
    # v_max           = 10815.

    # if H_mean = 18000000 (18e6)
    PACS_min = 30
    PACS_max = 365.0
    dt_orbit_min = 26500.0
    dt_orbit_max = 51000.0
    dt_eclipse_min = 2650.0
    dt_eclipse_max = 3750.0
    theta_slew_min = 0.00000875
    theta_slew_max = 0.00002100
    v_min = 3350.0
    v_max = 4950.0

    # Setup Var_Info

    var_info = {}

    # Save variable names in var_info
    var_info.update({"rand_inputs": rand_inputs})
    var_info.update({"fixed_inputs": fixed_inputs})
    var_info.update({"output_vars": output_vars})

    var_info.update(
        {
            # Save fixed variables
            "RE": RE,
            "mu": mu,
            "phi": phi,
            "c": c,
            "A_s": A_s,
            "sun_i": sun_i,
            "dt_slew": dt_slew,
            "M": M,
            "rho": rho,
            "A": A,
            "n": n,
            "omega_max": omega_max,
            "P_hold": P_hold,
            "I_d": I_d,
            "nu": nu,
            "LT": LT,
            "eps_deg": eps_deg,
            "r_lw": r_lw,
            "n_sa": n_sa,
            "rho_sa": rho_sa,
            "t": t,
            "D": D,
            "I_bodyx": I_bodyx,
            "I_bodyy": I_bodyy,
            "I_bodyz": I_bodyz,
            "I_max": I_max,
            "I_min": I_min,
            "theta": theta,
            "RD": RD,
            # Save mean and variance of random variables
            "H_mean": H_mean,
            "phi_mean": phi_mean,
            "Po_mean": Po_mean,
            "F_s_mean": F_s_mean,
            # 'theta_mean' : theta_mean,
            "L_sp_mean": L_sp_mean,
            "q_mean": q_mean,
            # 'RD_mean' : RD_mean,
            "L_a_mean": L_a_mean,
            "C_d_mean": C_d_mean,
            "H_std": H_std,
            "phi_std": phi_std,
            "Po_std": Po_std,
            "F_s_std": F_s_std,
            # 'theta_std' : theta_std,
            "L_sp_std": L_sp_std,
            "q_std": q_std,
            # 'RD_std' : RD_std,
            "L_a_std": L_a_std,
            "C_d_std": C_d_std,
            # Save min/max of coupling variables
            "PACS_min": PACS_min,
            "PACS_max": PACS_max,
            "dt_orbit_min": dt_orbit_min,
            "dt_orbit_max": dt_orbit_max,
            "dt_eclipse_min": dt_eclipse_min,
            "dt_eclipse_max": dt_eclipse_max,
            "theta_slew_min": theta_slew_min,
            "theta_slew_max": theta_slew_max,
            "v_min": v_min,
            "v_max": v_max,
        }
    )
    return var_info


def run(x, **kwargs):
    """Run full coupled fire satellite multiphysics problem

    Parameters
    ----------
    x : np.ndarray (6 or 7, n)
        Input design vars
        x[0] = H
        x[1] = phi
        x[2] = Po
        x[3] = F_s
        x[4] = L_sp
        x[5] = q
        x[6] = L_a
        x[7] = C_d

    **kwargs :
        Optional keyword arguments
        usehifi = (True, FALSE) select hifidelity calculations
        var_info = firesat.setup(), dictionary containing fixed parameters
        feedforward = (TRUE, False) select feedforward or feedback
            implementation

    Returns
    -------
    q : np.ndarray (3, n)
        Quantities of interest for system
        q[0] = P_tot
        q[1] = A_sa
        q[2] = tau_tot
    """
    # Get optional arguments
    usehifi = kwargs.get("hifi", False)  # use hifi calculations
    sat_params = kwargs.get("var_info", setup())  # fixed parameters
    feedforward = kwargs.get("feedforward", True)
    debug = kwargs.get("debug")

    # Compute orbit discipline
    x_orb = x[[0, 1]]
    q_orb = orbit(x_orb, sat_params)

    if feedforward:
        # Feed-forward implementation of Fire Satellite problem

        # Compute attitude control discipline
        x_atd = np.vstack((x[[0, 3, 4, 5, 6, 7]]))  # H, F_s, L_sp, q, L_a, C_d
        y_atd = np.vstack((q_orb[[0, 1, 3]]))  # v, dt_orbit, theta_slew
        q_atd = attitude(x_atd, y_atd, sat_params)

        # Compute power discipline
        x_pow = x[[2, 3]]  # Po, F_s
        y_pow = np.vstack((q_atd[1], q_orb[[1, 2]]))  # PACS, dt_orbit, dt_eclipse
        q_pow = power(x_pow, y_pow, sat_params, usehifi=usehifi)

    else:
        # Feedback coupling not implemented yet
        raise NotImplementedError

    if debug:
        print("\n")
        cplvars = {
            "PACS": q_atd[1],
            "dt_orbit": q_orb[1],
            "dt_eclipse": q_orb[2],
            "theta_slew": q_orb[3],
            "v": q_orb[0],
        }
        for var, vals in cplvars.items():
            print(f"{var:12s} min={vals.min():17.10f}  max={vals.max():17.10f}")

    n = x.shape[1]  # number of samples
    q = np.zeros((3, n))
    q[0] = q_pow[0]
    q[1] = q_pow[1]
    q[2] = q_atd[0]
    return q


if __name__ == "__main__":
    """
    Run full fire detection satellite coupled system

    Design variables:
        x = [H, Po, F_s, theta, L_sp, q, RD, L_a, C_d]

    Quantities of interest:
        q = [P_tot, A_sa, tau_tot]

    Command line parameters:

    $ python firesat.py n m

    with
        n = number of random samples (default n=1e6)
        m = random seed for reproducibility (default m=None)
    """

    # import sys
    # from timeit import default_timer as tic

    # # Set problem parameters
    # if len(sys.argv) >= 3:
    #     # Set random seed
    #     seed = int(float(sys.argv[2]))
    #     np.random.seed(seed)
    # if len(sys.argv) >= 2:
    #     # Set number of random samples
    #     n = int(float(sys.argv[1]))
    # else:
    #     n = int(1e4)

    # usehifi = True

    # # Setup satellite parameters
    # sat_params = setup()

    # # Sample iid input variables with mean and covariance matrices
    # input_vars = sat_params['rand_inputs']
    # x = mvn(input_vars, n)

    # # Run feedforward system
    # print(f'Run system with {n:,d} samples...', end=' ')
    # t = tic()
    # qoi = system(x, var_info=sat_params, feedforward=True, hifi=usehifi)
    # tf = tic()
    # print(f'{tf-t:.3f} sec')
    # print(qoi)

    # # # plot results
    # # import matplotlib.pyplot as plt
    # # import seaborn as sns
    # # sns.set(style="white", palette="muted", color_codes=True)
    # # qoi_opts = {}
    # # qoi_opts['latex'] = [r'$P_{tot}$', r'$A_{sa}$', r'$\tau_{tot}$']
    # # qoi_opts['color'] = ['b', 'r', 'g']

    # # f, axes = plt.subplots(1, 3) #, figsize=(7, 7), sharex=True)
    # # for i, qoivar in enumerate(range(3)):
    # #     sns.distplot(qoi[i], color=qoi_opts['color'][i], ax=axes[i])
    # #     axes[i].set_xlabel(qoi_opts['latex'][i])

    # # plt.setp(axes, yticks=[])
    # # plt.show()
