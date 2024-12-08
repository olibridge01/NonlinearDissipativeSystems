import random
import numpy as np

class Constants(object):
    """Physical constants for the MD/RPMD simulations."""

    angstrom = 1e-10
    hbar = 1.054571817e-34
    q_e = 1.602176634e-19
    m_e = 9.1093837015e-31
    a_0 = 5.291772109e-11
    kB = 3.166814535e-6
    E_h = ((hbar ** 2) / (m_e * (a_0 ** 2)))

    m_p = 1836.152672
    m_sys = m_p
    m_alpha = m_p

    w_b = (2 * np.pi * 1.49896229e+13) / (E_h / hbar)
    w_c = (2 * np.pi * 1.49896229e+13) / (E_h / hbar)
    V_0 = 4.14173925e-20 / E_h
    pot_1 = m_sys * (w_b ** 2)
    pot_2 = ((m_sys ** 2) * (w_b ** 4)) / (4 * V_0)
    min_val = -np.sqrt((4 * V_0) / (m_sys * (w_b ** 2)))


def inv_T(T: float) -> float:
    """Computes the inverse temperature (1/(kB*T))."""
    return 1 / (Constants.kB * T)


def init_p(p_dist: float, N: int, m: float, beta: float) -> float:
    """Initialise momenta from Gaussian distribution (used for periodic resampling of momenta from Boltzmann distribution)."""
    sigma_p = np.sqrt((N * m) / beta)
    p_dist = random.gauss(0, sigma_p)
    return p_dist


def rpmd_C(N: int) -> np.ndarray:
    """Generates the transformation matrix to transform the positions and momenta into the normal mode representation."""
    C = np.zeros((N, N))
    for j in range(N):
        for k in range(N):
            if k == 0:
                C[j][k] = np.sqrt(1 / N)
            elif 0 < k <= N / 2 - 1:
                C[j][k] = np.sqrt(2 / N) * np.cos((2 * np.pi * (j+1) * k) / N)
                if np.isclose(C[j][k],0):
                    C[j][k] = 0
            elif k == N / 2:
                C[j][k] = np.sqrt(1 / N) * np.power(-1, (j+1))
            elif k >= N / 2 + 1:
                C[j][k] = np.sqrt(2 / N) * np.sin((2 * np.pi * (j+1) * k) / N)
                if np.isclose(C[j][k],0):
                    C[j][k] = 0
    return C


def rpmd_E(omega_k: float, m: float, dt: float, constraint: bool = False) -> np.ndarray:
    """
    Takes value for omega_k and returns the corresponding evolution matrix for the symplectic integration.

    Args:
        omega_k (float): Normal mode frequency.
        m (float): Particle mass (atomic units).
        dt (float): Timestep (atomic units).
        constraint (bool, optional): Boolean to indicate whether system centroid is constrained/unconstrained. Defaults to False.
    """
    E = np.zeros((2, 2))
    if omega_k == 0:
        if constraint:
            E = np.identity(2) # Centroid mode is not evolved if constrained
        else:
            E[0][0] = 1
            E[0][1] = 0
            E[1][0] = dt / m
            E[1][1] = 1
    else:
        E[0][0] = np.cos(omega_k * dt)
        E[0][1] = -m * omega_k * np.sin(omega_k * dt)
        E[1][0] = (1 / (m * omega_k)) * np.sin(omega_k * dt)
        E[1][1] = np.cos(omega_k * dt)
    return E