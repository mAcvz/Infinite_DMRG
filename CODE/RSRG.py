from numba import jit 
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as plt
import numpy as np
import os
from numpy import linalg as LA
from module_work import generate_hamiltonian

# Number of iterations and spins in the system
n_iterations = 2
n_spins = 3
proj_states = 2**n_spins  # Number of projected states

# Initialize variables
n_jj = n_spins
delta = 1e-5
sigma_x = np.array([[0,1],[1,0]])
h_grid = np.linspace(-3, 0 , 50)
vec_gs = np.zeros(len(h_grid), dtype=float)

# Loop over the external field values
for ii, hh in enumerate(h_grid):
    # Generate the Hamiltonian matrix for the given number of spins and external field value
    Hamiltonian = generate_hamiltonian(n_spins, hh)
    
    # Get the eigenvalues and eigenstates of the Hamiltonian
    eigenvals, eigenstates = LA.eigh(Hamiltonian)
    
    # Initial steps for system doubling
    dim_H = Hamiltonian.shape[0]
    id_N_spins = np.eye(dim_H)
    id_N_spins_1 = np.eye(int(dim_H / 2))
    A_single_spin = np.kron(id_N_spins_1, sigma_x)
    B_single_spin = np.kron(sigma_x, id_N_spins_1)
    left_sys = np.kron(Hamiltonian, id_N_spins) 
    right_sys = np.kron(id_N_spins, Hamiltonian)
    interaction = np.kron(A_single_spin, B_single_spin)
    
    # Iterate through system doubling process
    n_jj = n_spins
    for jj in range(0, n_iterations):
        # Double the size of the system
        n_jj = n_jj * 2
        
        # Construct the new Hamiltonian
        Hamiltonian = left_sys + right_sys + interaction 
        eigenvals, eigenstates = LA.eigh(Hamiltonian)
        
        # Project onto a reduced space
        Pr = eigenstates[:, :proj_states]
        Pr_T = np.transpose(Pr)
        Pr_H = np.conj(Pr_T)
        Hamiltonian = Pr_H @ (Hamiltonian @ Pr)
        
        # Update system components based on the reduced space
        dim_H = Hamiltonian.shape[0]
        id_N_spins = np.eye(dim_H)
        A_single_spin = Pr_H @ (np.kron(id_N_spins, A_single_spin) @ Pr)
        B_single_spin = Pr_H @ (np.kron(B_single_spin, id_N_spins) @ Pr)
        left_sys = np.kron(Hamiltonian, id_N_spins) 
        right_sys = np.kron(id_N_spins, Hamiltonian)
        interaction = np.kron(A_single_spin, B_single_spin)
        prev_step_eng = eigenvals[0] / n_jj
    
    # Calculate the ground state energy for the current external field value
    eigenvals, eigenstates = LA.eigh(Hamiltonian)
    vec_gs[ii] = eigenvals[0]

# Save the ground state energies to a file
np.savetxt("GS_RSRG.txt", vec_gs)
