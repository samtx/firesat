import numpy as np
from firesat import system


def mvn(var_list, n=100000):
    """Generate iid multivariate normal random samples for variables in
    var_list, samples are generated independently with a diagonal covariance
    matrix

    Args:
        var_list : list of strings of variable names,
            e.g. var_list = ['H', 'phi'] for orbit
                var_list = ['H', 'F_s', 'L_sp', 'q', 'L_a', 'C_d'] for attitude
                var_list = ['Po', 'F_s'] for power
        n : int, number of random samples
    Returns:
        nd.array((d, n)) : numpy array of n samples for each d variable,
                        where d = len(var_list)
    """
    sat_params = system.setup()
    n_vars = len(var_list)
    MU = np.empty(n_vars)
    COV = np.zeros((n_vars, n_vars))
    for i, var in enumerate(var_list):
        MU[i] = sat_params[var + "_mean"]
        COV[i, i] = sat_params[var + "_std"] ** 2
    return np.random.multivariate_normal(MU, COV, int(n)).T


def uniform(var_list, n=100000):
    """Generate iid multivariate uniform random samples for variables in
    var_list, samples are generated independently.

    Args:
        var_list : list of strings of variable names,
                    e.g. var_list = ['H', 'phi'] for orbit
                        var_list = ['H', 'F_s', 'L_sp', 'q', 'L_a', 'C_d'], attitude
                        var_list = ['Po', 'F_s'] for power
        n : int, number of random samples
    Returns:
        nd.array((d, n)) : numpy array of n samples for each d variable, where
                           d = len(var_list)
    """
    sat_params = system.setup()
    n_vars = len(var_list)
    out = np.empty((n_vars, n))
    for i, var in enumerate(var_list):
        out[i] = np.random.uniform(
            sat_params[var + "_min"], sat_params[var + "_max"], n
        )
    return out
