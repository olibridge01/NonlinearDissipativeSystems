[![Python](https://img.shields.io/badge/Python-3776AB?logo=python&logoColor=fff)](#)
[![NumPy](https://img.shields.io/badge/NumPy-4DABCF?logo=numpy&logoColor=fff)](#)


# Quantum Rates in Dissipative Systems with Spatially Varying Friction

This repository contains a compact Python implementation of the [Ring Polymer Molecular Dynamics (RPMD)](https://doi.org/10.1063/1.1777575) method for computing approximate quantum rate constants in dissipative systems with linear and nonlinear friction. The code is used to simulate the dynamics of a ring polymer in a 1D potential, under the influence of a bath of harmonic oscillators to emulate a condensed phase environment.

RPMD results produced from this code contributed to the paper ["*Quantum Rates in Dissipative Systems with Spatially Varying Friction*", Bridge et al., J. Chem. Phys. (2024).](https://doi.org/10.1063/5.0216823)



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
- **`acf_example.ipynb`**: Example implementation of `autocorrelation.py`.


## Brief Background

<p align="center">
  <img src="https://github.com/olibridge01/NonlinearDissipativeSystems/assets/86416298/01d37871-359f-4006-a4b3-0a2783bd4e5b"  width="300">
</p>

The ring polymer Hamiltonian is given by
$$H_N(\mathbf{p},\mathbf{q}) = \sum_{j=1}^N \left[\frac{p_j^2}{2m} + \frac{1}{2}m\omega_N^2(q_j-q_{j+1})^2 + V(q_j)\right]$$

The RPMD flux-side time-correlation function is given by
$$\widetilde{C}_{fs}^{(N)}(t) = \frac{1}{(2\pi\hbar)^{Nf}} \int \text{d}^{Nf}\mathbf{p} \int \text{d}^{Nf}\mathbf{q} \ \text{e}^{-\beta_N H_N(\mathbf{p},\mathbf{q})} \delta[\bar{s}(\mathbf{q})] \bar{v}_s(\mathbf{p},\mathbf{q}) h[\bar{s}(\mathbf{q}_t)]$$



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