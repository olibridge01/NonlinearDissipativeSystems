import random
import numpy as np
from typing import List, Union, Tuple

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

def friction(eta_0: float, q: Union[np.ndarray, float]) -> Union[np.ndarray, float]:
    """
    Friction function for nonlinear bath coupling studied in this project.

    Args:
        eta_0 (float): Function parameter.
        q (float): System position.
    """
    return np.sqrt(eta_0) * q * (1 - np.exp(-(q ** 2) / 2))


def friction_deriv(eta_0: float, q: Union[np.ndarray, float]) -> Union[np.ndarray, float]:
    """Derivative of friction function for nonlinear bath coupling studied in this project."""
    return np.sqrt(eta_0) * (1 - np.exp(-(q ** 2) / 2) + (q ** 2) * np.exp(-(q ** 2) / 2))

def nm_transform(p: np.ndarray, x: np.ndarray, px: np.ndarray, RPMD_C: np.ndarray) -> np.ndarray:
    """Transform the positions and momenta into the normal mode representation."""
    px[:, 0] = np.dot(p, RPMD_C)
    px[:, 1] = np.dot(x, RPMD_C)
    return px

def rev_transform(px: np.ndarray, RPMD_C: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Transform the positions and momenta back into the original representation."""
    p = np.dot(RPMD_C, px[:, 0])
    x = np.dot(RPMD_C, px[:, 1])
    return p, x

def evolve(px: np.ndarray, RPMD_E: List[np.ndarray], RPMD_E_CONST: List[np.ndarray], constraint: bool = False) -> np.ndarray:
    """Evolve the positions and momenta in the velocity Verlet algorithm using the evolution matrix."""
    if constraint:
        for i_bead, row in enumerate(px):
            px[i_bead, :] = np.dot(RPMD_E_CONST[i_bead], row)
    else:
        for i_bead, row in enumerate(px):
            px[i_bead, :] = np.dot(RPMD_E[i_bead], row)

    return px

def heaviside(x: float) -> int:
    """Heaviside step function."""
    if x > 0:
        return 1
    elif x <= 0:
        return 0

def centroid(arr: Union[np.ndarray, float], n_beads: int) -> Union[np.ndarray, float]:
    """Compute the centroid of an array of bead positions."""
    if isinstance(arr, float):
        return arr
    centroid = 0
    for entry in arr:
        centroid += (1 / n_beads) * entry
    return centroid