import random
import numpy as np
from typing import List, Tuple, Union

from NonlinearDissipativeSystems.src.utils import *
from NonlinearDissipativeSystems.src.system import System, LinearSystem, NonlinearSystem
from NonlinearDissipativeSystems.src.bath import BathMode, LinearBathMode, NonlinearBathMode

def k_cl_tst(T: float) -> float:
    """Computes the rate constant using classical transition state theory."""
    return ((np.sqrt(2) * Constants.w_b) / (2 * np.pi)) * np.exp(-inv_T(T) * Constants.V_0)

def create_systembath(beta: float, n_bathmodes: int, gamma: float, n_beads: int, dt: float, linear: bool = True) -> Tuple[System, List[BathMode]]:
    """
    Initialise system and bath objects.
    
    Args:
        beta (float): Inverse temperature (1/(kB*T)).
        n_bathmodes (int): Number of bathmodes.
        gamma (float): Bath friction coefficient.
        n_beads (int): Number of ring polymer beads.
        dt (float): Timestep (atomic units).
        linear (bool): Boolean to indicate whether system-bath is linear or nonlinear.
    """
    # Instantiate array of BathMode instances
    bath = []
    for i_bathmode in range(n_bathmodes):
        if linear:
            bath.append(LinearBathMode(beta, n_bathmodes, i_bathmode, gamma, n_beads, dt))
        else:
            bath.append(NonlinearBathMode(beta, n_bathmodes, i_bathmode, gamma, n_beads, dt))

    # Instantiate system
    if linear:
        syst = LinearSystem(beta, n_beads, dt)
    else:
        syst = NonlinearSystem(beta, n_beads, dt)

    return syst, bath

def systembath_thermostat(syst: System, bath: List[BathMode], n_beads: int, beta: float) -> Tuple[System, List[BathMode]]:
    """Resample momenta from Boltzmann distribution."""
    syst.p = np.vectorize(init_p)(syst.p, n_beads, Constants.m_sys, beta)
    for bathmode in bath:
        bathmode.p = np.vectorize(init_p)(bathmode.p, n_beads, Constants.m_alpha, beta)
    return syst, bath


class RateCalc(object):
    """
    Rate calculation class that stores the parameters for a given calculation.
    """
    def __init__(
        self, 
        T: float, 
        n_bathmodes: int, 
        n_samp_kappa: int, 
        n_samp_fe: int, 
        n_samp_rd: int, 
        n_equil: int, 
        n_evol_kappa: int, 
        n_evol_fe: int,
        linear: bool = True
    ):
        """
        Args:
            T (float): Temperature (K).
            n_bathmodes (int): Number of bathmodes.
            n_samp_kappa (int): Number of samples for calculation of transmission coefficient (kappa).
            n_samp_fe (int): Number of samples for calculation of free energy.
            n_samp_rd (int): Number of samples for calculation of the reactant distribution.
            n_equil (int): Number of cycles in the equilibration phase (discarded before sampling phase).
            n_evol_kappa (int): Number of timesteps in the transmission evolution phase.
            n_evol_fe (int): Number of timesteps in the free energy evolution phase.
        
        Returns:
            None.
        """
        self.T = T
        self.beta = inv_T(self.T)
        self.n_bathmodes = n_bathmodes
        self.n_samp_kappa = n_samp_kappa
        self.n_samp_fe = n_samp_fe
        self.n_samp_rd = n_samp_rd
        self.n_equil = n_equil
        self.n_evol_kappa = n_evol_kappa
        self.n_evol_fe = n_evol_fe
        self.dt = 0.05 * ((2 * np.pi) / (-Constants.w_c * np.log(0.5 / self.n_bathmodes)))
        self.n_config = 50 # Number of Verlet iterations for constrained centroid sampling
        self.linear = linear

    def classical_transmission(self, gamma_factor: float, i_repeat: int) -> np.ndarray:
        """
        Produce the classical transmission coefficient for a given friction constant.
        
        Args:
            gamma_factor (float): Friction constant expressed in units of w_b.
            i_repeat (int): Repeat index.
        """
        gamma = gamma_factor * Constants.w_b

        # Define numerator and denominator of kappa as per Bennett-Chandler method
        numer = np.zeros(self.n_evol_kappa)
        denom = 0

        # Set n_beads to 1 (classical particle)
        n_beads = 1

        # Instantiate array of BathMode instances
        syst, bath = create_systembath(self.beta, self.n_bathmodes, gamma, n_beads, self.dt, linear=self.linear)
        print(f'{gamma} , no. {i_repeat} started.')

        # Loop over n_samples
        for i_sample in range(self.n_samp_kappa):

            # Initial positions/momenta
            syst.p = init_p(syst.p, n_beads, Constants.m_sys, self.beta)
            syst.x = 0
            for bathmode in bath:
                bathmode.p = init_p(bathmode.p, n_beads, Constants.m_alpha, self.beta)
                bathmode.classical_constrained_config()

            denom += (syst.p / Constants.m_sys) * heaviside(syst.p)
            print(f"Sample: {i_sample + 1}/{self.n_samp_kappa}   ", end='\r')
            A = syst.p / Constants.m_sys

            # Velocity verlet algorithm
            for i_evolution in range(self.n_evol_kappa):
                syst.p += -(self.dt / 2) * syst.sys_force(bath, gamma)
                for bathmode in bath:
                    bathmode.p += -(self.dt / 2) * bathmode.bathmode_force(syst.x)

                syst.x += self.dt * (syst.p / Constants.m_sys)
                for bathmode in bath:
                    bathmode.x += self.dt * (bathmode.p / Constants.m_alpha)

                syst.p += -(self.dt / 2) * syst.sys_force(bath, gamma)
                for bathmode in bath:
                    bathmode.p += -(self.dt / 2) * bathmode.bathmode_force(syst.x)

                B = heaviside(syst.x)
                numer[i_evolution] += A * B

        numer /= self.n_samp_kappa
        denom /= self.n_samp_kappa
        kappa = numer / denom
        return kappa

    def classical_mean_force(self, gamma: float, constrain_x: float, verbose: bool = False) -> float:
        """
        Compute the centroid mean force for a classical particle (N=1).

        Args:
            gamma (float): Friction constant.
            constrain_x (float): Position value of constrained particle.
        """
        n_beads = 1
        syst, bath = create_systembath(self.beta, self.n_bathmodes, gamma, n_beads, self.dt, linear=self.linear)
        syst.x = constrain_x
        mean_force = 0

        # Equilibration phase
        for i_equilibrium in range(self.n_equil):
            if verbose:
                print(f"Equil: {i_equilibrium + 1}/{self.n_equil}", end='\r')

            for bathmode in bath:
                bathmode.p = init_p(bathmode.p, n_beads, Constants.m_alpha, self.beta)

            # Velocity verlet algorithm
            for i_segment in range(self.n_evol_fe):
                for bathmode in bath:
                    bathmode.p += -(self.dt / 2) * bathmode.bathmode_constrained_force(constrain_x)

                for bathmode in bath:
                    bathmode.x += self.dt * (bathmode.p / Constants.m_alpha)

                for bathmode in bath:
                    bathmode.p += -(self.dt / 2) * bathmode.bathmode_constrained_force(constrain_x)

        # Sampling phase
        for i_sample in range(self.n_samp_fe):
            if verbose: 
                print(f"Sample: {i_sample + 1}/{self.n_samp_fe}", end='\r')

            for bathmode in bath:
                bathmode.p = init_p(bathmode.p, n_beads, Constants.m_alpha, self.beta)

            # Velocity verlet algorithm
            for i_evolution in range(self.n_evol_fe):
                for bathmode in bath:
                    bathmode.p += -(self.dt / 2) * bathmode.bathmode_constrained_force(constrain_x)

                for bathmode in bath:
                    bathmode.x += self.dt * (bathmode.p / Constants.m_alpha)

                for bathmode in bath:
                    bathmode.p += -(self.dt / 2) * bathmode.bathmode_constrained_force(constrain_x)

            mean_force += syst.sys_force(bath, gamma)

        mean_force /= self.n_samp_fe
        if verbose:
            print(f'x={constrain_x} done.')
        return mean_force

    def propagate_trajectory(self, syst: System, bath: List[BathMode], gamma: float, constraint: bool = False) -> Tuple[System, List[BathMode]]:
        """Propagate the system RP through one timestep."""

        syst.p += -(self.dt / 2) * syst.sys_force(bath, gamma)
        for bathmode in bath:
            bathmode.p += -(self.dt / 2) * bathmode.bathmode_force(syst.x)

        syst.px = nm_transform(syst.p, syst.x, syst.px, syst.RPMD_C)
        for bathmode in bath:
            bathmode.px = nm_transform(bathmode.p, bathmode.x, bathmode.px, bathmode.RPMD_C)

        syst.px = evolve(syst.px, syst.RPMD_E, syst.RPMD_E_CONST, constraint=constraint)
        for bathmode in bath:
            bathmode.px = evolve(bathmode.px, bathmode.RPMD_E, bathmode.RPMD_E_CONST)

        syst.p, syst.x = rev_transform(syst.px, syst.RPMD_C)
        for bathmode in bath:
            bathmode.p, bathmode.x = rev_transform(bathmode.px, bathmode.RPMD_C)

        syst.p += -(self.dt / 2) * syst.sys_force(bath, gamma)
        for bathmode in bath:
            bathmode.p += -(self.dt / 2) * bathmode.bathmode_force(syst.x)

        return syst, bath

    def rpmd_constrained_centroid(self, gamma: float, constrain_x: float, n_beads: int) -> Tuple[np.ndarray, np.ndarray]:
        """Compute a list of initial positions for the system and bath drawn from the constrained-centroid distribution."""
        
        sys_config = np.zeros((n_beads, self.n_samp_kappa))
        bath_config = np.zeros((n_beads, self.n_bathmodes, self.n_samp_kappa))

        syst, bath = create_systembath(self.beta, self.n_bathmodes, gamma, n_beads, self.dt, linear=self.linear)
        for i_bead in range(n_beads):
            syst.x[i_bead] = constrain_x

        for i_equilibrium in range(self.n_equil):
            # print(f'Equil: {i_equilibrium} / {n_equilibrium}')
            syst, bath = systembath_thermostat(syst, bath, n_beads, self.beta)
            for i_segment in range(self.n_evol_kappa):
                syst, bath = self.propagate_trajectory(syst, bath, gamma, constraint=True)

        for i_sample in range(self.n_samp_kappa):
            syst, bath = systembath_thermostat(syst, bath, n_beads, self.beta)
            for i_step_config in range(self.n_config):
                syst, bath = self.propagate_trajectory(syst, bath, gamma, constraint=True)

            sys_config[:, i_sample] = syst.x
            for i_bathmode, bathmode in enumerate(bath):
                bath_config[:, i_bathmode, i_sample] = bathmode.x

        return sys_config, bath_config

    def rpmd_transmission(self, gamma_factor: float, n_beads: int, i_repeat: int, verbose: bool = False) -> np.ndarray:
        """Compute RPMD transmission factor."""

        gamma = gamma_factor * Constants.w_b
        numer = np.zeros(self.n_evol_kappa)
        denom = 0
        syst, bath = create_systembath(self.beta, self.n_bathmodes, gamma, n_beads, self.dt, linear=self.linear)

        if verbose:
            print(f'{gamma}, no. {i_repeat} started.')

        sys_config, bath_config = self.rpmd_constrained_centroid(gamma, 0, n_beads)

        for i_sample in range(self.n_samp_kappa):
            if verbose:
                print(f'Sample: {i_sample} / {self.n_samp_kappa}', end='\r')

            syst.p = np.vectorize(init_p)(syst.p, n_beads, Constants.m_sys, self.beta)
            syst.x = sys_config[:, i_sample]

            for i_bathmode, bathmode in enumerate(bath):
                bathmode.p = np.vectorize(init_p)(bathmode.p, n_beads, Constants.m_sys, self.beta)
                bathmode.x = bath_config[:, i_bathmode, i_sample]

            p_centroid = centroid(syst.p, n_beads)

            denom += (p_centroid / Constants.m_sys) * heaviside(p_centroid)
            A = p_centroid / Constants.m_sys

            for i_evolution in range(self.n_evol_kappa):
                syst, bath = self.propagate_trajectory(syst, bath, gamma)

                x_centroid = centroid(syst.x, n_beads)
                B = heaviside(x_centroid)
                numer[i_evolution] += A * B

        numer /= self.n_samp_kappa
        denom /= self.n_samp_kappa
        kappa = numer / denom
        return kappa

    def rpmd_mean_force(self, gamma: float, n_beads: int, constrain_x: float, verbose: bool = False) -> float:
        """Compute RPMD centroid-constrained mean force."""
        mean_force = 0
        syst, bath = create_systembath(self.beta, self.n_bathmodes, gamma, n_beads, self.dt, linear=self.linear)
        for i_bead in range(n_beads):
            syst.x[i_bead] = constrain_x

        for i_equilibrium in range(self.n_equil):
            if verbose:
                print(f'Equil: {i_equilibrium} / {self.n_equil}', end='\r')

            syst, bath = systembath_thermostat(syst, bath, n_beads, self.beta)
            for i_segment in range(self.n_evol_fe):
                syst, bath = self.propagate_trajectory(syst, bath, gamma, constraint=True)

        for i_sample in range(self.n_samp_fe):
            if verbose:
                print(f'Sample: {i_sample} / {self.n_samp_fe}', end='\r')
            syst, bath = systembath_thermostat(syst, bath, n_beads, self.beta)
            for i_evolution in range(self.n_evol_fe):
                syst, bath = self.propagate_trajectory(syst, bath, gamma, constraint=True)

            for i_bead in range(n_beads):
                mean_force += (1 / n_beads) * syst.sys_force(bath, gamma)[i_bead]

        mean_force /= self.n_samp_fe
        return mean_force

    def rpmd_mf_array(self, gamma_factor: float, n_beads: int, n_points: int, i_repeat: int, min: float = -2.0, max: float = 0.0) -> np.ndarray:
        """
        Compute series of mean forces along the reaction coordinate.
        
        Args:
            gamma_factor (float): Friction constant expressed in units of w_b.
            n_beads (int): Number of RP beads.
            n_points (int): Number of points along the reaction coordinate.
            i_repeat (int): Repeat index.
            min (float, optional): Minimum value of reaction coordinate.
            max (float, optional): Maximum value of reaction coordinate.
        """
        gamma = gamma_factor * Constants.w_b
        x = np.linspace(min, max, n_points) # Reaction coordinate values
        mf = np.zeros(n_points)

        for idx, value in enumerate(x):
            mf[idx] = self.rpmd_mean_force(gamma=gamma, n_beads=n_beads, constrain_x=value)
            print(f'Repeat {i_repeat}: {value} done.')

        return mf

    def reactant_distribution(self, gamma_factor: float, n_beads: int, i_repeat: int, verbose: bool = True) -> np.ndarray:
        """Compute reactant distribution of RP."""

        if verbose:
            print(f'{gamma_factor}, no. {i_repeat} started.')

        gamma = gamma_factor * Constants.w_b
        centroid_positions = np.zeros(self.n_samp_rd)
        syst, bath = create_systembath(self.beta, self.n_bathmodes, gamma, n_beads, self.dt, linear=self.linear)

        for i_bead in range(n_beads):
            syst.x[i_bead] = Constants.min_val

        for i_equilibrium in range(self.n_equil):
            if verbose:
                print(f'Equil: {i_equilibrium} / {self.n_equil}', end='\r')

            syst, bath = systembath_thermostat(syst, bath, n_beads, self.beta)
            for i_segment in range(self.n_evol_fe):
                syst, bath = self.propagate_trajectory(syst, bath, gamma)

        for i_sample in range(self.n_samp_rd):
            if verbose:
                print(f'Sample: {i_sample} / {self.n_samp_rd}', end='\r')

            syst, bath = systembath_thermostat(syst, bath, n_beads, self.beta)
            for i_evolution in range(self.n_evol_fe):
                syst, bath = self.propagate_trajectory(syst, bath, gamma)

            if centroid(syst.x, n_beads) > 0:
                syst.x = -syst.x
                for bathmode in bath:
                    bathmode.x = -bathmode.x

            x_centroid = centroid(syst.x, n_beads)
            centroid_positions[i_sample] += x_centroid

        return centroid_positions