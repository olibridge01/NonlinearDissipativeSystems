![Python3.10](https://img.shields.io/badge/python-3.8%20%7C%203.9%20%7C%203.10-blue)

![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)
![NumPy](https://img.shields.io/badge/numpy-%23013243.svg?style=for-the-badge&logo=numpy&logoColor=white)

# Quantum Reaction Rates in Nonlinear Dissipative Systems
## Part III Project

**Oli Bridge** (<ob344@cam.ac.uk>) - *St Catharine's College, Cambridge*

Part III Project in the *Althorpe group* - http://www-stuart.ch.cam.ac.uk/index.html

---

<img src="https://github.com/olibridge01/NonlinearDissipativeSystems/assets/86416298/7b7e55d3-136a-4483-aa91-531cd155e891"  width="400">

### Abstract

Low energy excitations of electron-hole pairs in metals induce spatially-dependent non-adiabatic effects (NAEs) by coupling the nuclear and electronic degrees of freedom. It has recently been proposed that one can model this with a system nonlinearly coupled to a bath of harmonic oscillators. Additionally, systems involving light atoms exhibit significant nuclear quantum effects (NQEs), which may be accounted for with ring polymer molecular dynamics (RPMD), a method for simulating approximate quantum dynamics that conserves the quantum Boltzmann distribution. In this work, we investigate how the quantum rate in a one-dimensional double-well potential is affected by position-dependent friction. Numerical calculations show that the recrossing dynamics, which behave classically at room temperature, are significantly affected by the nonlinear friction. For strong coupling, quantum effects on the rate are also observed since the free energy contribution from the bath becomes strongly dependent on position. In the deep-tunnelling regime, a systematic lowering of the instanton crossover temperature with bath friction is observed to heavily impact the free energy barriers for both linear and nonlinear coupling. Extension of this model to a one-dimensional treatment of hydrogen hopping in bulk palladium shows that the inclusion of NAEs yields increases in the rate, primarily due to altered recrossing, that are greater than from the inclusion of NQEs.

---

## Setup and Installation
Running the following command will install the required Python packages for executing this code:
```
python setup.py
```

## Directory Structure

```
.
├── engine/
│   ├── autocorrelation.py
│   ├── bath.py
│   └── utils.py
├── figures/
├── scripts/
│   ├── example1/
│   └── example2/
├── thesis/
├── acf_example.ipynb
└── README.md
```

- **`engine/`**: Contains the main code for the RPMD simulations.

  - **`autocorrelation.py`**: Code for computing position autocorrelation functions of a ring polymer in a 1D potential.
  - **`bath.py`**: Code for linear/nonlinear bath rate calculations.
  - **`utils.py`**: General functions for constant temperature RPMD simulations.
 
- **`figures/`**: Thesis figures.

- **`scripts/`**: Scripts for running simulations.

- **`thesis/`**: Full project thesis and slides for presentation of the project.

- **`acf_example.ipynb`**: Example implementation of `autocorrelation.py`.

---
## Brief Background

The ring polymer Hamiltonian is given by
$$H_N(\textbf{p},\textbf{q}) = \sum_{j=1}^N \left[\frac{p_j^2}{2m} + \frac{1}{2}m\omega_N^2(q_j-q_{j+1})^2 + V(q_j)\right]$$

The RPMD flux-side time-correlation function is given by
$$\widetilde{C}_{fs}^{(N)}(t) = \frac{1}{(2\pi\hbar)^{Nf}} \int \text{d}^{Nf}\textbf{p} \int \text{d}^{Nf}\textbf{q} \ \text{e}^{-\beta_N H_N(\textbf{p},\textbf{q})} \delta[\bar{s}(\textbf{q})] \bar{v}_s(\textbf{p},\textbf{q}) h[\bar{s}(\textbf{q}_t)]$$



