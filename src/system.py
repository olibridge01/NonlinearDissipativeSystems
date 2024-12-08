import random
import numpy as np
from typing import List

from NonlinearDissipativeSystems.src.utils import Constants, friction_deriv, rpmd_C, rpmd_E

class System(object):
    """
    Base system object.
    """
    def __init__(self, beta: float, n: int, dt: float):
        """
        Args:
            beta (float): Inverse temperature (1/(kB*T)).
            n (int): Number of ring polymer beads.
            dt (float): Timestep (atomic units).
        """
        # Inverse temperature
        self.beta = beta

        # Initialise system position/momentum
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

    def sys_force(self, bath: List, gamma: float) -> np.ndarray:
        """Computes the force on the system DoF from the potential and the bath."""
        raise NotImplementedError


class LinearSystem(System):
    """
    Linear system object that stores its position/momentum, as well as a force method.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def sys_force(self, bath: List, gamma: float) -> np.ndarray:
        """Computes the force on the system DoF from the potential and the bath."""

        # Interchange force if using different potentials
        force = -(Constants.pot_1 * self.x) + (Constants.pot_2 * (self.x ** 3)) + ((2 * Constants.m_sys * gamma *
                                                                                    Constants.w_c) / np.pi) * self.x
        
        # Add force from each bath mode
        for bathmode in bath:
            force -= bathmode.g_alpha * bathmode.x
        return force
    

class NonlinearSystem(System):
    """
    Nonlinear system object that stores its position/momentum, as well as a force method.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def sys_force(self, bath: List, gamma: float) -> np.ndarray:
        """Computes the force on the system DoF from the potential and the nonlinear bath."""
        force = -(Constants.pot_1 * self.x) + (Constants.pot_2 * (self.x ** 3)) + ((2 * Constants.w_c) / np.pi) * friction_deriv(Constants.m_sys*gamma, self.x)
        for bathmode in bath:
            force -= bathmode.g_alpha * friction_deriv(Constants.m_sys * gamma, self.x) * bathmode.x
        return force