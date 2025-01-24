- This package allows the simulation of the scattering intensity and form factor of particles in the shape of a right prism with a regular polygon cross section with any number of sides. Both 1D intensity profile and 2D images can be simulated.
- Parameters for describing the shape are the length $L$ of the right prism, the number $n$ of the regular polygon cross section and the edge length $E$ of the polygon. Alternatively, the average radius $R_{ave}$ can be used instead of $E$.
- It can also be used with the ipywidgets package to manually fit experimental data with the provided models.

### Installation
Package installation can be done using pip command: 
```
pip install prismformfactors
```
as the package is available on PyPI (https://pypi.org/project/prismformfactors).

### Example
add .dat file 

and jupyter notebook (with many comments)

make the modification about $\sigma_L$ and $\sigma_w$ in vizualize


### References
- Jules Marcone, Jaime Gabriel Trazo, Rahul Nag, Claire Goldmann, Nicolas
Ratel-Ramond, Cyrille Hamon and Marianne Imperor-Clerc, "Form factor of prismatic particles for Small Angle Scattering analysis", Journal of Applied Crystallography (2025)

- Theory behind code:
    Wuttke "Numerically stable form factor of any polygon and polyhedron." Journal of Applied Crystallography 54.2 (2021): 580-587
    https://doi.org/10.1107/S1600576721001710

Created on September 5th, 2024 and revised on January 24th, 2025

Creator: Marianne Impéror and Jules Marcone

Original Contributors: Marianne Impéror and Jules Marcone
