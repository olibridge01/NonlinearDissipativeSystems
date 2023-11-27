# Part III Project - Quantum Reaction Rates in Nonlinear Dissipative Systems
---

Oli Bridge, <ob344@cam.ac.uk>

St Catharine's College, Cambridge

Part III Project in Althorpe group.

---

## Setup and Installation
Running the following command will install the required Python packages for executing this code:
```
python setup.py
```

## Directories

- `engine`: Main programs that run the simulations
- `tools`: Secondary functions
- `scripts`: Plotters for the data produced by my code etc.
- `data`: Results and plots from path integral simulations

The ring polymer Hamiltonian is given by
$$H_N(\textbf{p},\textbf{q}) = \sum_{j=1}^N \left[\frac{p_j^2}{2m} + \frac{1}{2}m\omega_N^2(q_j-q_{j+1})^2 + V(q_j)\right]$$

The RPMD flux-side time-correlation function is given by
$$\widetilde{C}_{fs}^{(N)}(t) = \frac{1}{(2\pi\hbar)^{Nf}} \int \text{d}^{Nf}\textbf{p} \int \text{d}^{Nf}\textbf{q} \ \text{e}^{-\beta_N H_N(\textbf{p},\textbf{q})} \delta[\bar{s}(\textbf{q})] \bar{v}_s(\textbf{p},\textbf{q}) h[\bar{s}(\textbf{q}_t)]$$
