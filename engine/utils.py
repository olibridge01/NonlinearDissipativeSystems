import numpy as np
import random

class Constants(object):
    """
    Physical constants for the MD/RPMD simulations.
    """
    # Atomic units
    angstrom = 1e-10
    hbar = 1.054571817e-34
    q_e = 1.602176634e-19
    m_e = 9.1093837015e-31
    a_0 = 5.291772109e-11
    kB = 3.166814535e-6
    E_h = ((hbar ** 2) / (m_e * (a_0 ** 2)))

    # Masses, temperature
    m_p = 1836.152672
    m_sys = m_p
    m_alpha = m_p

    # External potential constants
    w_b = (2 * np.pi * 1.49896229e+13) / (E_h / hbar)
    w_c = (2 * np.pi * 1.49896229e+13) / (E_h / hbar)
    V_0 = 4.14173925e-20 / E_h
    pot_1 = m_sys * (w_b ** 2)
    pot_2 = ((m_sys ** 2) * (w_b ** 4)) / (4 * V_0)
    min_val = -np.sqrt((4 * V_0) / (m_sys * (w_b ** 2)))


def inv_T(T):
    """
    Computes the inverse temperature (1/(kB*T))

    Parameters
    ----------
    T : float
        Temperature (K).
    
    Returns
    -------
    inv_T : float
        Inverse temperature (1/(kB*T)) in atomic units.
    """
    return 1 / (Constants.kB * T)


def init_p(p_distribution, n, m, beta):
    """
    Initialise momenta from gaussian distribution (used for periodic resampling of momenta from Boltzmann distribution).
    
    Parameters
    ----------
    p_distribution : float
        Initial momenta value.
    n : int
        Number of ring polymer beads.
    m : float
        Particle mass in atomic units.
    beta : float
        Inverse temperature (1/(kB*T)).
    
    Returns
    -------
    p_distribution : float
        Updated momenta value.
    """
    sigma_p = np.sqrt((n * m) / beta)
    p_distribution = random.gauss(0, sigma_p)
    return p_distribution


def rpmd_C(n):
    """
    Generates the transformation matrix to transform the positions and momenta into the normal mode representation.

    Parameters
    ----------
    n : int
        Number of ring polymer beads.
    
    Returns
    -------
    C : array
        Transformation matrix with shape (N,N).
    """
    C = np.zeros((n, n))
    for j in range(n):
        for k in range(n):
            if k == 0:
                C[j][k] = np.sqrt(1 / n)
            elif 0 < k <= n / 2 - 1:
                C[j][k] = np.sqrt(2 / n) * np.cos((2 * np.pi * (j+1) * k) / n)
                if np.isclose(C[j][k],0):
                    C[j][k] = 0
            elif k == n / 2:
                C[j][k] = np.sqrt(1 / n) * np.power(-1, (j+1))
            elif k >= n / 2 + 1:
                C[j][k] = np.sqrt(2 / n) * np.sin((2 * np.pi * (j+1) * k) / n)
                if np.isclose(C[j][k],0):
                    C[j][k] = 0
    return C


def rpmd_E(omega_k, mass, dt, constraint=False):
    """
    Takes value for omega_k and returns the corresponding evolution matrix for the symplectic integration.

    Parameters
    ----------
    omega_k : float
        Normal mode frequency.
    dt : float
        Timestep (atomic units).
    constraint : bool
        Boolean to indicate whether system centroid is constrained or unconstrained.

    Returns
    -------
    E : array
        Evolution matrix with shape (2,2).
    """
    E = np.zeros((2, 2))
    if omega_k == 0:
        # Deals with the omega = 0 limit of the evolution matrix, which would be undefined using numpy functions
        if constraint:
            E = np.identity(2)
        else:
            E[0][0] = 1
            E[0][1] = 0
            E[1][0] = dt / mass
            E[1][1] = 1
    else:
        E[0][0] = np.cos(omega_k * dt)
        E[0][1] = -mass * omega_k * np.sin(omega_k * dt)
        E[1][0] = (1 / (mass * omega_k)) * np.sin(omega_k * dt)
        E[1][1] = np.cos(omega_k * dt)
    return E