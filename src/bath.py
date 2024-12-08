import random
import numpy as np
from typing import List, Tuple, Union

from NonlinearDissipativeSystems.src.utils import Constants, inv_T, init_p, rpmd_C, rpmd_E, friction

class BathMode(object):
    """
    Base bathmode object.
    """
    def __init__(self, beta: float, n_bathmodes: int, i_bathmode: int, gamma: float, n: int, dt: float):
        """
        Args:
            beta (float): Inverse temperature (1/(kB*T)).
            n_bathmodes (int): Number of bathmodes.
            i_bathmode (int): Index of bathmode.
            gamma (float): Bath friction coefficient.
            n (int): Number of ring polymer beads.
            dt (float): Timestep (atomic units).
        """
        self.beta = beta
        self.bathmode_id = i_bathmode
        self.gamma = gamma

        # Set bathmode frequency and coupling parameter
        self.w_alpha = -Constants.w_c * np.log((self.bathmode_id + 0.5) / n_bathmodes)
        self.g_alpha = self.w_alpha * np.sqrt((2 * self.gamma * Constants.m_sys * Constants.m_alpha * Constants.w_c) /
                                              (np.pi * n_bathmodes))
        
        # Initialise bathmode position/momentum
        self.p = np.zeros(n)
        self.x = np.zeros(n)
        self.px = np.zeros((n, 2))

        # Initialise RPMD transformation and evolution matrices for velocity verlet algorithm
        self.RPMD_C = rpmd_C(n)
        self.omega_k = np.zeros(n)
        for i_bead in range(n):
            self.omega_k[i_bead] = 2 * (n / self.beta) * np.sin((i_bead * np.pi) / n)
        self.RPMD_E = [rpmd_E(self.omega_k[i_bead], Constants.m_sys, dt) for i_bead in range(n)]
        self.RPMD_E_CONST = [rpmd_E(self.omega_k[i_bead], Constants.m_sys, dt, constraint=True) for i_bead in range(n)]

    def bathmode_force(self, sys_x: np.ndarray) -> np.ndarray:
        """Computes the force on the bathmode DoF from the potential and the system."""
        raise NotImplementedError

    def bathmode_constrained_force(self, constrain_x: float) -> float:
        """Compute the force on the bathmode DoF with a constraint on the system position."""
        raise NotImplementedError


class LinearBathMode(BathMode):
    """
    Bathmode object that stores its position/momentum/frequency etc.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def bathmode_force(self, sys_x: np.ndarray) -> np.ndarray:
        """Computes the force on the bathmode DoF from the potential and the system."""

        force = Constants.m_alpha * (self.w_alpha ** 2) * self.x - self.g_alpha * sys_x
        return force

    def classical_constrained_config(self):
        """Sample initial bath mode position from harmonic Gaussian distribution (FOR CLASSICAL ONLY)."""

        self.x = np.random.normal(0, np.sqrt(1 / (self.beta * Constants.m_alpha * (self.w_alpha ** 2))))
        return

    def bathmode_constrained_force(self, constrain_x: float) -> float:
        """Compute the force on the bathmode DoF with a constraint on the system position (FOR CLASSICAL ONLY)."""

        force = Constants.m_alpha * (self.w_alpha ** 2) * self.x - self.g_alpha * constrain_x
        return force


class NonlinearBathMode(LinearBathMode):
    """
    Bath mode object that stores its position/momentum/frequency etc.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def bathmode_force(self, sys_x: np.ndarray) -> np.ndarray:
        """Computes the force on the bathmode DoF from the potential and the system with nonlinear coupling."""

        force = Constants.m_alpha * (self.w_alpha ** 2) * self.x - self.g_alpha * friction(Constants.m_sys * self.gamma, sys_x)
        return force

    def bathmode_constrained_force(self, constrain_x):
        """Compute the force on the bathmode DoF with a constraint on the system position (FOR CLASSICAL ONLY)."""

        force = Constants.m_alpha * (self.w_alpha ** 2) * self.x - self.g_alpha * friction(Constants.m_sys * self.gamma, constrain_x)
        return force