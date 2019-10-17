import numpy as np


def power(x, y, var_info, fidelity=0, **kwargs):
    """Power model to compute the size of solar array required to power
    satellite. Takes the power requirements from the attitude model and the
    orbit and eclipse time from the orbit model. Outputs the moments of inertia
    based on the computed solar array geometry.

    Parameters
    ----------
    x : np.ndarray (2, n)
        Input design variables
        x[0] = Po
        x[1] = F_s

    y : np.ndarray (3, n)
        Input coupling variables
        y[0] = PACS
        y[1] = dt_orbit
        y[2] = dt_eclipse

    var_info : dict
        Dictionary containing fixed parameters for problem

    **kwargs :
        Optionary keyword arguments to select high fidelity computation
        fidelity : (int)
            0 = lowest fidelity  (default)
            1 = medium fidelity
            2 = high fidelity

    Returns
    -------
    q : np.ndarray (4, n)
        Quantities of interest and output coupling variables
        q[0] = P_tot
        q[1] = A_sa
        q[2] = I_max
        q[3] = I_min
    """
    # Setup fixed input variables
    sun_i = var_info["sun_i"]
    I_d = var_info["I_d"]
    nu = var_info["nu"]
    LT = var_info["LT"]
    eps_deg = var_info["eps_deg"]
    r_lw = var_info["r_lw"]
    n_sa = var_info["n_sa"]
    rho_sa = var_info["rho_sa"]
    t = var_info["t"]
    D = var_info["D"]
    I_bodyx = var_info["I_bodyx"]
    I_bodyy = var_info["I_bodyy"]
    I_bodyz = var_info["I_bodyz"]

    # Unpack design variables
    Po, F_s = x

    # Unpack coupling variables
    PACS, dt_orbit, dt_eclipse = y

    # Compute output quantitites
    P_tot = PACS + Po
    T_e = dt_eclipse
    T_d = dt_orbit - T_e
    P_e = P_tot
    P_d = P_tot
    P_sa = (P_e * T_e / 0.6 + P_d * T_d / 0.8) / T_d
    PBOL = nu * F_s * I_d * np.cos(sun_i)
    PEOL = PBOL * (1 - eps_deg) ** (LT)
    A_sa = P_sa / PEOL

    # Compute solar array geometry quantities
    L = np.sqrt(A_sa * r_lw / n_sa)
    W = np.sqrt(A_sa / (r_lw * n_sa))

    if fidelity > 0:
        # use high fidelity computation
        out = inertia_high_fidelity(W, L, D, t, rho_sa, n_sa)
        I_sax, I_say, I_saz, m_sa = out

    else:
        # use low fidelity calculation
        out = inertia_low_fidelity(W, L, D, t, rho_sa, n_sa)
        I_sax, I_say, I_saz, m_sa = out

    # print(f'I_sa x, y, z = {I_sax.mean()}, {I_say.mean()}, {I_saz.mean()}')

    I_sax += I_bodyx
    I_say += I_bodyy
    I_saz += I_bodyz

    I_tot = np.array([I_sax, I_say, I_saz]).T
    if len(I_tot.shape) > 1:
        I_max = np.amax(I_tot, axis=1)
        I_min = np.amin(I_tot, axis=1)
    else:
        I_max = np.amax(I_tot, axis=0)
        I_min = np.amin(I_tot, axis=0)

    # Assemble Power outputs
    n = x.shape[1]  # number of samples
    q = np.zeros((4, n))
    q[0] = P_tot
    q[1] = A_sa
    q[2] = I_max
    q[3] = I_min
    return q


def inertia_low_fidelity(W, L, D, t, rho_sa, n_sa):
    # ** double check the math here... there might be more solar panels than are assumed.
    m_sa = rho_sa * L * W * t
    I_sax = m_sa * (1 / 12 * (L ** 2 + t ** 2) + (D + L / 2) ** 2)
    I_say = m_sa / 12 * (t ** 2 + W ** 2)
    I_saz = m_sa * (1 / 12 * (L ** 2 + W ** 2) + (D + L / 2) ** 2)
    return I_sax, I_say, I_saz, m_sa


def inertia_high_fidelity(W, L, D, t, rho_sa, n_sa):
    """Hi-fidelity inertia tensor calculation using tetgen to
    create a 3D tetrahedron mesh of the solar panels. Then
    compute the intertia tensor from the sum of the element
    interias
    """
    import tetgen

    I_sax = np.zeros(W.size)
    I_say = np.zeros(W.size)
    I_saz = np.zeros(W.size)
    m_sa = np.zeros(W.size)
    for i in range(W.size):
        Wi = W[i]
        Li = L[i]
        # set vertices
        verts = np.array(
            [
                [-Wi / 2, D, t / 2],
                [Wi / 2, D, t / 2],
                [Wi / 2, D, -t / 2],
                [-Wi / 2, D, -t / 2],
                [-Wi / 2, D + Li, t / 2],
                [Wi / 2, D + Li, t / 2],
                [Wi / 2, D + Li, -t / 2],
                [-Wi / 2, D + Li, -t / 2],
            ]
        )
        # set facets
        facets = np.array(
            [
                [0, 1, 2],
                [0, 2, 3],
                [0, 1, 5],
                [0, 5, 4],
                [0, 3, 7],
                [0, 4, 7],
                [1, 2, 5],
                [2, 5, 6],
                [3, 2, 6],
                [3, 6, 7],
                [4, 5, 6],
                [4, 6, 7],
            ]
        )
        tet = tetgen.TetGen(verts, facets)
        nodes, elems = tet.tetrahedralize(quality=1, verbose=0)
        # Compute inertia tensor from tetrahedrons
        nelem = elems.shape[0]
        dIxx = np.zeros(nelem)
        dIyy = np.zeros(nelem)
        dIzz = np.zeros(nelem)
        dMass = np.zeros(nelem)
        for j, el in enumerate(elems):
            # Gather element nodes
            pts = nodes[el]
            a = pts[0]
            b = pts[1]
            c = pts[2]
            d = pts[3]
            # Find centroid x, y, z coordinates
            xx = (a[0] + b[0] + c[0] + d[0]) / 4.0
            yy = (a[1] + b[1] + c[1] + d[1]) / 4.0
            zz = (a[2] + b[2] + c[2] + d[2]) / 4.0
            # Find element volume
            dV = (
                1
                / 6.0
                * abs(
                    (a[0] - d[0])
                    * ((b[1] - d[1]) * (c[2] - d[2]) - (b[2] - d[2]) * (c[1] - d[1]))
                    + (a[1] - d[1])
                    * ((b[2] - d[2]) * (c[0] - d[0]) - (b[0] - d[0]) * (c[2] - d[2]))
                    + (a[2] - d[2])
                    * ((b[0] - d[0]) * (c[1] - d[1]) - (b[1] - d[1]) * (c[0] - d[0]))
                )
            )
            # Compute element inertia
            dMass[j] = rho_sa * dV
            dIxx[j] = rho_sa * dV * (yy * yy + zz * zz)
            dIyy[j] = rho_sa * dV * (xx * xx + zz * zz)
            dIzz[j] = rho_sa * dV * (xx * xx + yy * yy)
        I_sax[i] = np.sum(dIxx) ** 2
        I_say[i] = np.sum(dIyy) ** 2
        I_saz[i] = np.sum(dIzz) ** 2
        m_sa[i] = np.sum(dMass)

    return I_sax, I_say, I_saz, m_sa
