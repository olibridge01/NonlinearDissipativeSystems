[![Python](https://img.shields.io/badge/Python-3776AB?logo=python&logoColor=fff)](#)
[![NumPy](https://img.shields.io/badge/NumPy-4DABCF?logo=numpy&logoColor=fff)](#)


# Quantum Rates in Dissipative Systems with Spatially Varying Friction

This repository contains a compact Python implementation of the [Ring Polymer Molecular Dynamics (RPMD)](https://doi.org/10.1063/1.1777575) method for computing approximate quantum rate constants in dissipative systems with linear and nonlinear friction. The code is used to simulate the dynamics of a ring polymer in a 1D potential, under the influence of a bath of harmonic oscillators to emulate a condensed phase environment.

RPMD results produced from this code contributed to the paper ["*Quantum Rates in Dissipative Systems with Spatially Varying Friction*", Bridge et al., J. Chem. Phys. (2024)](https://doi.org/10.1063/5.0216823). For experimental data (including ML-MCTDH results), and plotting scripts for the figures in the paper, please see the paper's [GitLab repository](https://gitlab.com/litman90/quantum_spatially_varying_friction).


## Installation
To clone the repository, run the following:
  ```bash
  git clone https://github.com/olibridge01/NonlinearDissipativeSystems.git
  cd NonlinearDissipativeSystems
  ```
For package management, set up a conda environment (called `nl_rpmd` as an example) and install the required packages as follows:
  ```bash
  conda create -n nl_rpmd python=3.10 anaconda
  conda activate nl_rpmd
  pip install -r requirements.txt
  ```

## Running the Code
The code is designed to be run from the command line. To run a simulation, use the `run_experiment.py` script in the `scripts/` directory. For example, to run a simulation with a linear friction bath, use the following command:
  ```bash
  python scripts/run_experiment.py -t [temp] --n_bathmodes [nmodes] -g [gamma] -n [nbeads] -c config.yaml -l
  ```
where `[temp]` is the temperature (K), `[nmodes]` is the number of bath modes, `[gamma]` is the friction coefficient, `[nbeads]` is the number of ring polymer beads, `config.yaml` is the configuration file, and `-l` is a flag that specifies linear friction. 
For a nonlinear bath, omit the `-l` flag:
  ```bash
  python scripts/run_experiment.py -t [temp] --n_bathmodes [nmodes] -g [gamma] -n [nbeads] -c config.yaml
  ```


## Directory Structure

The directory structure is as follows:
```
.
├── figs/
├── results/
├── src/
│   ├── autocorrelation.py
│   ├── bath.py
│   ├── rate.py
│   ├── system.py
│   └── utils.py
├── scripts/
│   ├── run_experiment.py
│   └── plotters.py
├── thesis/
├── acf_example.ipynb
└── README.md
```

- **`figs/`**: Thesis figures.
- **`results/`**: Stores simulation results.
- **`src/`**: RPMD code.
  - `autocorrelation.py`: Code for reproducing the [original RPMD paper](https://doi.org/10.1063/1.1777575) results.
  - `bath.py`: Linear and nonlinear `BathMode()` classes.
  - `rate.py`: RPMD rate constant calculation.
  - `system.py`: Linear and nonlinear `System()` classes.
  - `utils.py`: General functions for constant temperature RPMD simulations.

- **`scripts/`**: Scripts for running simulations, plotting etc.
- **`thesis/`**: Contains my Part III thesis for reference.
- **`acf_example.ipynb`**: Example implementation of `autocorrelation.py`.


## Brief Background

### Ring Polymer Hamiltonian

<p align="center">
  <img src="https://github.com/olibridge01/NonlinearDissipativeSystems/assets/86416298/01d37871-359f-4006-a4b3-0a2783bd4e5b"  width="300">
</p>

The ring polymer Hamiltonian for a 1D potential is given by

$$
H_N(\mathbf{p},\mathbf{q}) = \sum_{j=1}^N \left[\frac{p_j^2}{2m} + \frac{1}{2}m\omega_N^2(q_j-q_{j+1})^2 + V(q_j)\right]
$$
where $N$ is the number of beads, $m$ is the mass of each bead, $\omega_N=N/\beta \hbar$, $\beta=1/k_BT$, and $p_j$ and $q_j$ are the momentum and position of the $j^{th}$ bead.

### RPMD Rates
The $N$-bead ring polymer approximation to the quantum rate constant is given by

$$
    k^{(N)}(T) = \frac{1}{Q_r^{(N)}(T)}\lim_{t\to\infty}\widetilde{C}^{(N)}_{fs}(t),
$$

where $\widetilde{C}_{fs}^{(N)}(t)$ is the ring polymer flux-side time-correlation function (TCF), and $Q_r^{(N)}(T)$ is the ring polymer partition function per unit volume.

The [Bennett-Chandler method](https://doi.org/10.1063/1.2883593) is a useful method for computing rate constants when the event of a reactive crossing is rare. We start by expressing the RPMD rate coefficient in terms of thermal averages:

$$
    k^{(N)}(T) = \lim_{t\to\infty}\frac{\langle \delta(q_s^{\ddagger} - \bar{q}_s)(\bar{p}_s/m) h(\bar{q}_s(t) - q_s^{\ddagger}) \rangle}{\langle h(q_s^{\ddagger} - \bar{q}_s)\rangle}.
$$

This can be written as a product of two terms:

$$
    k^{(N)}(T) = \kappa(\tau_p) k^{\text{QTST}}(T),
$$

where $\tau_p$ is the plateau time. The transmission coefficient, $\kappa(t)$, is defined as:

$$
    \kappa(t) = \frac{\langle \delta(q_s^{\ddagger} - \bar{q}_s)(\bar{p}_s/m_s) h(\bar{q}_s(t) - q_s^{\ddagger}) \rangle}{\langle \delta(q_s^{\ddagger} - \bar{q}_s) (\bar{p}_s/m_s) h(\bar{p}_s)\rangle},
$$

where we define $q_s$ ad the system reaction coordinate, $q_s^{\ddagger}$ as the transition state, $m_s$ as the system mass, and $\bar{q}_s$ and $\bar{p}_s$ as the centroid position and momentum of the ring polymer, respectively. The RPMD TST rate constant, $k^{\text{QTST}}(T)$, is given by

$$
    k^{\text{QTST}}(T) = \frac{1}{(2\pi \beta m_s)^{1/2}}\:p(q_s^0)\:\exp\left(-\beta\int_{q_s^0}^{q_s^{\ddagger}}\text{d}q_s'\:\frac{\text{d}\mathcal{F}(q_s')}{\text{d}q_s'}\right),
$$

where $\mathcal{F}(q_s)$ is the free energy along the reaction coordinate, and $p(q_s^0)$ is the probability of the ring polymer centroid being at $q_s^0$, a point in the reactant well. For more details on the rate calculation mathematical background, see my [Part III thesis](thesis/partIII_thesis.pdf).


## Citation

**Oli Bridge** (<olibridge@rocketmail.com>) - *St Catharine's College, Cambridge*

This work is from my Part III project in the [Althorpe group](http://www-stuart.ch.cam.ac.uk/index.html) at the [University of Cambridge](https://www.cam.ac.uk/), and has contributed to a [publication](https://doi.org/10.1063/5.0216823) in the Journal of Chemical Physics. If you use this code, please cite the following:

```
@article{bridge_quantum_2024,
	title = {Quantum rates in dissipative systems with spatially varying friction},
	volume = {161},
	issn = {0021-9606},
	url = {https://doi.org/10.1063/5.0216823},
	doi = {10.1063/5.0216823},
	pages = {024110},
	number = {2},
	journaltitle = {The Journal of Chemical Physics},
	author = {Bridge, Oliver and Lazzaroni, Paolo and Martinazzo, Rocco and Rossi, Mariana and Althorpe, Stuart C. and Litman, Yair},
	date = {2024-07},
}
```