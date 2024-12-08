import random
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple

from NonlinearDissipativeSystems.engine.utils import init_p, rpmd_C, rpmd_E


#--------------------Frequencies--------------------

def omegas(N: int, beta: float) -> np.ndarray:
    """Generates the ring polymer normal mode frequencies."""
    omega_N = N / beta
    omegas = np.zeros(N)
    for i in range(N):
        omegas[i] = 2 * omega_N * np.sin(i * np.pi / N)
    return omegas

#------------------Autocorrelation Function------------------

class AutoCorrelation(object):
    def __init__(
        self, 
        force: callable,
        beta: float = 1.0,
        mass: float = 1.0,
        dt: float = 0.05,
        n_samp: int = 1000,
        n_equil: int = 100,
        n_evol: int = 500
    ):
        """
        Args:
            force (function): Function for the force acting on the particle.
            beta (float): Inverse temperature (1/(kB*T)).
            mass (float): Particle mass in atomic units.
            dt (float): Timestep length in atomic units.
            n_samp (int): Number of samples to take in the sampling phase.
            n_equil (int): Number of cycles in the equilibration phase (discarded before sampling phase).
            n_evol (int): Number of timesteps in the evolution phase.
        
        Returns:
            None.
        """
        self.mass = mass
        self.beta = beta
        self.dt = dt
        self.n_samp = n_samp
        self.n_evol = n_evol
        self.n_equil = n_equil
        self.force = force

    def classical_verlet_step(self, x: float, p: float) -> tuple:
        """Velocity Verlet algorithm for a classical 1D particle."""

        p += (self.dt / 2) * self.force(x)
        x += self.dt * (p / self.mass)
        p += (self.dt / 2) * self.force(x)
        return x, p

    def classical_autocorrelation(self) -> np.ndarray:
        """Compute <x(0)x(t)> for a classical particle in a given 1D potential."""

        # Set number of beads to 1 (classical particle) and initialise array of x(0)x(t) values
        N = 1
        xx = np.zeros(self.n_evol)

        # Initialise particle momentum and position
        p = 0
        x = 0

        # Equilibriation phase
        for i_equilibrium in range(self.n_equil):
            p = init_p(p, N, self.mass, self.beta) # Momentum resampled every cycle to avoid nonergodicity

            # Velocity verlet algorithm
            for i_evol in range(self.n_evol):
                x, p = self.classical_verlet_step(x, p)

        # Sampling phase
        for i_sample in range(self.n_samp):
            # Resampling momentum each cycle
            p = init_p(p, N, self.mass, self.beta)

            # Setting A to current x value i.e. x(0)
            A = x

            # Velocity verlet
            for i_evol in range(self.n_evol):
                x, p = self.classical_verlet_step(x, p)

                # Setting B to current x value i.e. x(t)
                B = x

                # Adding to x(0)x(t) array for corresponding time value
                xx[i_evol] += A * B

        xx /= self.n_samp
        return xx

    def rpmd_verlet_step(self, N: int, x: np.ndarray, p: np.ndarray, RPMD_C: np.ndarray, omegas: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Velocity Verlet algorithm for a ring polymer particle in a 1D potential.
        
        Args:
            N (int): Number of ring polymer beads.
            x (array): Array of bead positions.
            p (array): Array of bead momenta.
            RPMD_C (array): Transformation matrix to transform the positions and momenta into the normal mode representation.
            omegas (array): Array of ring polymer normal mode frequencies.
        """
        p += (self.dt / 2) * self.force(x)

        px_vectors = np.zeros((N, 2))
        px_vectors[:, 0] = np.dot(p, RPMD_C)
        px_vectors[:, 1] = np.dot(x, RPMD_C)

        for i in range(N):
            px_vectors[i, :] = np.dot(rpmd_E(omegas[i], self.mass, self.dt), px_vectors[i, :])

        p = np.dot(RPMD_C, px_vectors[:, 0])
        x = np.dot(RPMD_C, px_vectors[:, 1])

        p += (self.dt / 2) * self.force(x)

        return x, p

    def rpmd_autocorrelation(self, N: int) -> np.ndarray:
        """Compute <x_N(0)x_N(t)> for an N-bead ring polymer in a given 1D potential."""

        # Create matrix objects for the trajectory propagator
        RPMD_C = rpmd_C(N)
        omega_list = omegas(N, self.beta)

        # Initialise array for the x(0)x(t) data, as well as the (p,x) vectors for each bead to be used in the normal mode basis
        xx = np.zeros(self.n_evol)
        px_vectors = np.zeros((N, 2))

        # Initialise arrays to store momentum and position of the beads
        p = np.zeros(N)
        x = np.zeros(N)

        # Initialise x
        for i in range(N):
            x[i] = 0

        # Equilibration phase
        for i_equilibrium in range(self.n_equil):

            # Resample momenta from Boltzmann distribution
            for i in range(N):
                p[i] = init_p(p[i], N, self.mass, self.beta)

            # Velocity Verlet with normal mode transformations
            for i_evol in range(self.n_evol):
                x, p = self.rpmd_verlet_step(N, x, p, RPMD_C, omega_list)

        # Sampling phase
        for i_sample in range(self.n_samp):

            # Resample momenta from Boltzmann distribution
            for i in range(N):
                p[i] = init_p(p[i], N, self.mass, self.beta)

            A = 0
            for i in range(N):
                A += x[i]
            A /= N

            # Velocity verlet
            for i_evol in range(self.n_evol):
                x, p = self.rpmd_verlet_step(N, x, p, RPMD_C, omega_list)

                B = 0
                for i in range(N):
                    B += x[i]
                B /= N

                # Add x(0)x(t) value to the xx array
                xx[i_evol] += (A * B)

        xx /= self.n_samp
        return xx