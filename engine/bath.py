import numpy as np
import random
from engine.utils import Constants, inv_T, init_p, rpmd_C, rpmd_E

def k_cl_tst(T):
    """
    Computes harmonic approximation to the rate constant using classical transition state theory.

    Parameters
    ----------
    T : float
        Temperature (K).

    Returns
    -------
    k_cl_tst : float
        Classical TST rate constant in atomic units.
    """
    return ((np.sqrt(2) * Constants.w_b) / (2 * np.pi)) * np.exp(-inv_T(T) * Constants.V_0)

#----------------------Classes----------------------

class LinearSystem:
    """
    Linear system object that stores its position/momentum, as well as a force method.
    """
    def __init__(self, beta, n, dt):
        """
        Parameters
        ----------
        beta : float
            Inverse temperature (1/(kB*T)).
        n : int
            Number of ring polymer beads.
        dt : float
            Timestep (atomic units).
        
        Returns
        -------
        None.
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

    def sys_force(self, bath, gamma):
        """
        Computes the force on the system DoF from the potential and the bath.

        Parameters
        ----------
        bath : array
            Array of BathMode objects.
        gamma : float
            Bath friction coefficient.

        Returns
        -------
        force : array
            Force on the system.
        """
        force = -(Constants.pot_1 * self.x) + (Constants.pot_2 * (self.x ** 3)) + ((2 * Constants.m_sys * gamma *
                                                                                    Constants.w_c) / np.pi) * self.x
        for bathmode in bath:
            force -= bathmode.g_alpha * bathmode.x
        return force
    

class NonlinearSystem(LinearSystem):
    """
    Nonlinear system object that stores its position/momentum, as well as a force method.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def sys_force(self, bath, gamma):
        """
        Computes the force on the system DoF from the potential and the nonlinear bath.

        Parameters
        ----------
        bath : array
            Array of BathMode objects. 
        gamma : float
            Bath friction coefficient.
        
        Returns
        -------
        force : array
            Force on the system.
        """
        force = -(Constants.pot_1 * self.x) + (Constants.pot_2 * (self.x ** 3)) + ((2 * Constants.w_c) / np.pi) * friction_deriv(Constants.m_sys*gamma, self.x)
        for bathmode in bath:
            force -= bathmode.g_alpha * friction_deriv(Constants.m_sys * gamma, self.x) * bathmode.x
        return force


class LinearBathMode:
    """
    Bathmode object that stores its position/momentum/frequency etc.
    """
    def __init__(self, beta, n_bathmodes, i_bathmode, gamma, n, dt):
        """
        Parameters
        ----------
        beta : float
            Inverse temperature (1/(kB*T)).
        n_bathmodes : int
            Number of bathmodes.
        i_bathmode : int
            Index of bathmode.
        gamma : float
            Bath friction coefficient.
        n : int
            Number of ring polymer beads.
        dt : float
            Timestep (atomic units).

        Returns
        -------
        None.
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

    def bathmode_force(self, sys_x):
        """
        Computes the force on the bathmode DoF from the potential and the system.
        
        Parameters
        ----------
        sys_x : array
            System position.
        
        Returns
        -------
        force : array
            Force on the bathmode.
        """
        force = Constants.m_alpha * (self.w_alpha ** 2) * self.x - self.g_alpha * sys_x
        return force

    def classical_constrained_config(self):
        """
        Sample initial bath mode position from harmonic gaussian distribution (FOR CLASSICAL ONLY)
        
        Returns
        -------
        None.
        """
        self.x = np.random.normal(0, np.sqrt(1 / (self.beta * Constants.m_alpha * (self.w_alpha ** 2))))
        return

    def bathmode_constrained_force(self, constrain_x):
        """
        Compute the force on the bathmode DoF with a constraint on the system position (FOR CLASSICAL ONLY).

        Parameters
        ----------
        constrain_x : float
            Position value of constrained system DoF.

        Returns
        -------
        force : float
            Force on the bathmode.
        """
        force = Constants.m_alpha * (self.w_alpha ** 2) * self.x - self.g_alpha * constrain_x
        return force


class NonlinearBathMode(LinearBathMode):
    """
    Bath mode object that stores its position/momentum/frequency etc.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def bathmode_force(self, sys_x):
        """
        Computes the force on the bathmode DoF from the potential and the system with nonlinear coupling.
        
        Parameters
        ----------
        sys_x : array
            System position.
        
        Returns
        -------
        force : array
            Force on the bathmode.
        """
        force = Constants.m_alpha * (self.w_alpha ** 2) * self.x - self.g_alpha * friction(Constants.m_sys * self.gamma, sys_x)
        return force

    def bathmode_constrained_force(self, constrain_x):
        """
        Compute the force on the bathmode DoF with a constraint on the system position (FOR CLASSICAL ONLY).

        Parameters
        ----------
        constrain_x : float
            Position value of constrained system DoF.

        Returns
        -------
        force : float
            Force on the bathmode.
        """
        force = Constants.m_alpha * (self.w_alpha ** 2) * self.x - self.g_alpha * friction(Constants.m_sys * self.gamma, constrain_x)
        return force

#---------------------Functions---------------------

def nm_transform(p, x, px, RPMD_C):
    """
    Transform the positions and momenta into the normal mode representation.
    
    Parameters
    ----------
    p : array
        Array of bead momenta.
    x : array
        Array of bead positions.
    px : array
        Array of bead positions/momenta in normal mode representation.
    RPMD_C : array
        Transformation matrix of shape (N,N) to transform positions/momenta into the normal mode basis.

    Returns
    -------
    px : array
        Updated array of bead positions/momenta in normal mode representation.
    """
    px[:, 0] = np.dot(p, RPMD_C)
    px[:, 1] = np.dot(x, RPMD_C)
    return px


def rev_transform(px, RPMD_C):
    """
    Transform the positions and momenta back into the original representation.
    
    Parameters
    ----------
    px : array
        Array of bead positions/momenta in normal mode representation.
    RPMD_C : array
        Transformation matrix of shape (N,N) to transform positions/momenta into the normal mode basis.
    
    Returns
    -------
    p : array
        Array of bead momenta.
    x : array
        Array of bead positions.
    """
    p = np.dot(RPMD_C, px[:, 0])
    x = np.dot(RPMD_C, px[:, 1])
    return p, x


def evolve(px, RPMD_E, RPMD_E_CONST, constraint=False):
    """
    Evolve the positions and momenta in the velocity Verlet algorithm using the evolution matrix.
    
    Parameters
    ----------
    px : array
        Array of bead positions/momenta in normal mode representation.
    RPMD_E : array
        Array of evolution matrices of shape (N,2,2) for each bead.
    RPMD_E_CONST : array
        Array of evolution matrices of shape (N,2,2) for each bead with the constraint applied.
    constraint : bool
        Boolean to indicate whether system centroid is constrained or unconstrained.

    Returns
    -------
    px : array
        Updated array of bead positions/momenta in normal mode representation.
    """
    if constraint:
        for i_bead, row in enumerate(px):
            px[i_bead, :] = np.dot(RPMD_E_CONST[i_bead], row)
    else:
        for i_bead, row in enumerate(px):
            px[i_bead, :] = np.dot(RPMD_E[i_bead], row)

    return px


def create_systembath(beta, n_bathmodes, gamma, n_beads, dt, linear=True):
    """
    Initialise system and bath objects.
    
    Parameters
    ----------
    beta : float
        Inverse temperature (1/(kB*T)).
    n_bathmodes : int
        Number of bathmodes.
    gamma : float
        Bath friction coefficient.
    n_beads : int
        Number of ring polymer beads.
    dt : float
        Timestep (atomic units).
    linear : bool
        Boolean to indicate whether system-bath is linear or nonlinear.

    Returns
    -------
    syst : System object
        System object.
    bath : array
        Array of BathMode objects.
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


def systembath_thermostat(syst, bath, n_beads, beta):
    """
    Resample momenta from Boltzmann distribution.

    Parameters
    ----------
    syst : System object
        System object.
    bath : array
        Array of BathMode objects.
    n_beads : int
        Number of ring polymer beads.
    beta : float
        Inverse temperature (1/(kB*T)).
    
    Returns
    -------
    syst : System object
        Updated system object.
    bath : array
        Updated array of BathMode objects.
    """
    syst.p = np.vectorize(init_p)(syst.p, n_beads, Constants.m_sys, beta)
    for bathmode in bath:
        bathmode.p = np.vectorize(init_p)(bathmode.p, n_beads, Constants.m_alpha, beta)
    return syst, bath


def heaviside(x):
    """
    Heaviside step function.
    """
    if x > 0:
        return 1
    elif x <= 0:
        return 0


def centroid(arr, n_beads):
    """
    Compute the centroid of an array of bead positions.
    
    Parameters
    ----------
    arr : array
        Array of bead positions.
    n_beads : int
        Number of ring polymer beads.
    
    Returns
    -------
    centroid : float
        Centroid position.
    """
    if isinstance(arr, float):
        return arr
    centroid = 0
    for entry in arr:
        centroid += (1/n_beads) * entry
    return centroid


def friction(eta_0, q):
    """
    Friction function for nonlinear bath coupling studied in this project.

    Parameters
    ----------
    eta_0 : float
        Function parameter.
    q : float
        System position.
    
    Returns
    -------
    friction : float
        Friction value.
    """
    return np.sqrt(eta_0) * q * (1 - np.exp(-(q ** 2) / 2))


def friction_deriv(eta_0, q):
    """
    Derivative of friction function for nonlinear bath coupling studied in this project.
    
    Parameters
    ----------
    eta_0 : float
        Function parameter.
    q : float
        System position.
    
    Returns
    -------
    friction_deriv : float
        Friction derivative.
    """
    return np.sqrt(eta_0) * (1 - np.exp(-(q ** 2) / 2) + (q ** 2) * np.exp(-(q ** 2) / 2))


class RateCalc:
    """
    Rate calculation object that stores the parameters for a given calculation.
    """
    def __init__(self, T, n_bathmodes, n_samp_kappa, n_samp_fe, n_samp_rd, n_equil, n_evol_kappa, n_evol_fe):
        """
        Parameters
        ----------
        T : float
            Temperature (K).
        n_bathmodes : int
            Number of bathmodes.
        n_samp_kappa : int
            Number of samples for calculation of transmission coefficient (kappa).
        n_samp_fe : int
            Number of samples for calculation of free energy.
        n_samp_rd : int
            Number of samples for calculation of the reactant distribution.
        n_equil : int
            Number of cycles in the equilibration phase (discarded before sampling phase).
        n_evol_kappa : int
            Number of timesteps in the transmission evolution phase.
        n_evol_fe : int
            Number of timesteps in the free energy evolution phase.
        
        Returns
        -------
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

    def classical_transmission(self, gamma_factor, i_repeat):
        """
        Produce the classical transmission coefficient for a given friction constant.
        
        Parameters
        ----------
        gamma_factor : float
            Friction constant expressed in units of w_b.
        i_repeat : int
            Repeat index.

        Returns
        -------
        kappa : array
            Array of transmission values corresponding to each timestep.
        """
        gamma = gamma_factor * Constants.w_b

        # Define numerator and denominator of kappa as per Bennett-Chandler method
        numer = np.zeros(self.n_evol_kappa)
        denom = 0

        # Set n_beads to 1 (classical particle)
        n_beads = 1

        # Instantiate array of BathMode instances
        syst, bath = create_systembath(self.beta, self.n_bathmodes, gamma, n_beads, self.dt)
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

    def classical_mean_force(self, gamma, constrain_x):
        """
        Compute the centroid mean force for a classical particle (N=1).
        
        Parameters
        ----------
        gamma : float
            Friction constant.
        constrain_x : float
            Position value of constrained particle.
        
        Returns
        -------
        mean_force : float
            Centroid mean force.
        """
        n_beads = 1
        syst, bath = create_systembath(self.beta, self.n_bathmodes, gamma, n_beads, self.dt)
        syst.x = constrain_x
        mean_force = 0

        for i_equilibrium in range(self.n_equil):
            # print(f"Equil: {i_equilibrium + 1}/{n_equilibrium}")

            for bathmode in bath:
                bathmode.p = init_p(bathmode.p, n_beads, Constants.m_alpha, self.beta)

            for i_segment in range(self.n_evol_fe):
                for bathmode in bath:
                    bathmode.p += -(self.dt / 2) * bathmode.bathmode_constrained_force(constrain_x)

                for bathmode in bath:
                    bathmode.x += self.dt * (bathmode.p / Constants.m_alpha)

                for bathmode in bath:
                    bathmode.p += -(self.dt / 2) * bathmode.bathmode_constrained_force(constrain_x)

        for i_sample in range(self.n_samp_fe):
            # print(f"Sample: {i_sample + 1}/{n_samples}")

            for bathmode in bath:
                bathmode.p = init_p(bathmode.p, n_beads, Constants.m_alpha, self.beta)

            for i_evolution in range(self.n_evol_fe):
                for bathmode in bath:
                    bathmode.p += -(self.dt / 2) * bathmode.bathmode_constrained_force(constrain_x)

                for bathmode in bath:
                    bathmode.x += self.dt * (bathmode.p / Constants.m_alpha)

                for bathmode in bath:
                    bathmode.p += -(self.dt / 2) * bathmode.bathmode_constrained_force(constrain_x)

            mean_force += syst.sys_force(bath, gamma)

        mean_force /= self.n_samp_fe
        print(f'x={constrain_x} done.')
        return mean_force

    def propagate_trajectory(self, syst, bath, gamma, constraint=False):
        """
        Propagate the system RP through one timestep.
        
        Parameters
        ----------
        syst : System object
            System object.
        bath : array
            Array of BathMode objects.
        gamma : float
            Bath friction coefficient.
        constraint : bool
            Boolean to indicate whether system centroid is constrained or unconstrained.
        
        Returns
        -------
        syst : System object
            Updated system object.
        """
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

    def rpmd_constrained_centroid(self, gamma, constrain_x, n_beads):
        """
        Compute a list of initial positions for the system and bath drawn from the constrained-centroid distribution.
        
        Parameters
        ----------
        gamma : float
            Friction constant.
        constrain_x : float
            Position value of constrained particle.
        n_beads : int
            Number of ring polymer beads.
        
        Returns
        -------
        sys_config : array
            Array of system positions.
        """
        n_config = 50

        sys_config = np.zeros((n_beads, self.n_samp_kappa))
        bath_config = np.zeros((n_beads, self.n_bathmodes, self.n_samp_kappa))

        syst, bath = create_systembath(self.beta, self.n_bathmodes, gamma, n_beads, self.dt)
        for i_bead in range(n_beads):
            syst.x[i_bead] = constrain_x

        for i_equilibrium in range(self.n_equil):
            # print(f'Equil: {i_equilibrium} / {n_equilibrium}')
            syst, bath = systembath_thermostat(syst, bath, n_beads, self.beta)
            for i_segment in range(self.n_evol_kappa):
                syst, bath = self.propagate_trajectory(syst, bath, gamma, constraint=True)

        for i_sample in range(self.n_samp_kappa):
            syst, bath = systembath_thermostat(syst, bath, n_beads, self.beta)
            for i_step_config in range(n_config):
                syst, bath = self.propagate_trajectory(syst, bath, gamma, constraint=True)

            sys_config[:, i_sample] = syst.x
            for i_bathmode, bathmode in enumerate(bath):
                bath_config[:, i_bathmode, i_sample] = bathmode.x

        return sys_config, bath_config

    def rpmd_transmission(self, gamma_factor, n_beads, i_repeat):
        """
        Compute RPMD transmission factor.
        
        Parameters
        ----------
        gamma_factor : float
            Friction constant expressed in units of w_b.
        n_beads : int
            Number of RP beads.
        i_repeat : int
            Repeat index.
        
        Returns
        -------
        kappa : array
            Array of transmission values corresponding to each timestep.
        """
        gamma = gamma_factor * Constants.w_b
        numer = np.zeros(self.n_evol_kappa)
        denom = 0
        syst, bath = create_systembath(self.beta, self.n_bathmodes, gamma, n_beads, self.dt)
        print(f'{gamma} , no. {i_repeat} started.')

        sys_config, bath_config = self.rpmd_constrained_centroid(gamma, 0, n_beads)

        for i_sample in range(self.n_samp_kappa):
            # print(f'Sample: {i_sample} / {n_samples}')

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

    def rpmd_mean_force(self, gamma, n_beads, constrain_x):
        """
        Compute RPMD centroid-constrained mean force.
        
        Parameters
        ----------  
        gamma : float
            Friction constant.
        n_beads : int
            Number of RP beads.
        constrain_x : float
            Position value of constrained particle.
        
        Returns
        -------
        mean_force : float
            Centroid mean force.
        """
        mean_force = 0
        syst, bath = create_systembath(self.beta, self.n_bathmodes, gamma, n_beads, self.dt)
        for i_bead in range(n_beads):
            syst.x[i_bead] = constrain_x

        for i_equilibrium in range(self.n_equil):
            # print(f'Equil: {i_equilibrium} / {n_equilibrium}')
            syst, bath = systembath_thermostat(syst, bath, n_beads, self.beta)
            for i_segment in range(self.n_evol_fe):
                syst, bath = self.propagate_trajectory(syst, bath, gamma, constraint=True)

        for i_sample in range(self.n_samp_fe):
            print(f'Sample: {i_sample} / {self.n_samp_fe}')
            syst, bath = systembath_thermostat(syst, bath, n_beads, self.beta)
            for i_evolution in range(self.n_evol_fe):
                syst, bath = self.propagate_trajectory(syst, bath, gamma, constraint=True)

            for i_bead in range(n_beads):
                mean_force += (1 / n_beads) * syst.sys_force(bath, gamma)[i_bead]

        mean_force /= self.n_samp_fe
        return mean_force

    def rpmd_mf_array(self, gamma_factor, n_beads, n_points, i_repeat, min=-2.0, max=0.0):
        """
        Compute series of mean forces along the reaction coordinate.
        
        Parameters
        ----------
        gamma_factor : float
            Friction constant expressed in units of w_b.
        n_beads : int
            Number of RP beads.
        n_points : int
            Number of points along the reaction coordinate.
        i_repeat : int
            Repeat index.
        min : float
            Minimum value of reaction coordinate.
        max : float
            Maximum value of reaction coordinate.
        
        Returns
        -------
        mf : array
            Array of mean force values corresponding to each reaction coordinate value. 
        """
        gamma = gamma_factor * Constants.w_b
        x = np.linspace(min, max, n_points)
        mf = np.zeros(n_points)

        for idx, value in enumerate(x):
            mf[idx] = self.rpmd_mean_force(gamma=gamma, n_beads=n_beads, constrain_x=value)
            print(f'Repeat {i_repeat}: {value} done.')

        return mf

    def reactant_distribution(self, gamma_factor, n_beads, i_repeat):
        """
        Compute reactant distribution of RP.
        
        Parameters
        ----------
        gamma_factor : float
            Friction constant expressed in units of w_b.
        n_beads : int
            Number of RP beads.
        i_repeat : int
            Repeat index.
        
        Returns
        -------
        centroid_positions : array
            Array of centroid positions sampled from the reactant distribution.
        """
        print(f'{gamma_factor} , no. {i_repeat} started.')

        gamma = gamma_factor * Constants.w_b
        centroid_positions = np.zeros(self.n_samp_rd)
        syst, bath = create_systembath(self.beta, self.n_bathmodes, gamma, n_beads, self.dt)

        for i_bead in range(n_beads):
            syst.x[i_bead] = Constants.min_val

        for i_equilibrium in range(self.n_equil):
            print(f'Equil: {i_equilibrium} / {self.n_equil}')

            syst, bath = systembath_thermostat(syst, bath, n_beads, self.beta)
            for i_segment in range(self.n_evol_fe):
                syst, bath = self.propagate_trajectory(syst, bath, gamma)

        for i_sample in range(self.n_samp_rd):
            print(f'Sample: {i_sample} / {self.n_samp_rd}')

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
