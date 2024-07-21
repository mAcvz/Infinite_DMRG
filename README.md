# Real Space and Density Matrix Renormalization Group Methods

This repository contains the implementation and results of Assignment 8 for the Quantum Information and Computing course at Università degli Studi di Padova. The focus is on the Real Space Renormalization Group (RSRG) and the Density Matrix Renormalization Group (DMRG) methods.

## Introduction
The Real Space Renormalization Group (RSRG) and Density Matrix Renormalization Group (DMRG) are powerful techniques in theoretical physics and quantum mechanics used to study complex systems. These methods offer crucial insights into fields such as condensed matter physics and quantum chemistry, enabling a deeper understanding of the behavior of matter at different scales.

## Real Space Renormalization Group (RSRG)

### Algorithm
The RSRG method approximates the ground state of a system by focusing on its low-energy states. The algorithm involves:
1. Creating the Hamiltonian for a system of size  $N$
2. Diagonalizing the Hamiltonian to find its eigenvalues and eigenvectors.
3. Constructing a projected Hamiltonian for the system of size $2N$ using the eigenstates.
4. Iterating these steps until the desired system size or convergence is reached.


## Density Matrix Renormalization Group (DMRG)

### Algorithm
DMRG is a modified RG algorithm with improved truncation rules for higher precision at the cost of slower growth in system size. The system size increases linearly instead of exponentially with each iteration. The algorithm involves:
1. Initializing operators for a system composed of four parts.
2. Enlarging the left and right blocks by adding sites.
3. Building and diagonalizing the Hamiltonian for 2m + 2  particles.
4. Constructing and projecting the density matrix for the ground state.


