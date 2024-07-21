import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as LA
import time
from initialize import initialize_operators
#
h_grid = np.linspace(-3, 0 , 50)
# Pauli matrices
sigma_x = np.array([[0,1],[1,0]])
sigma_z = np.array([[1,0],[0,-1]])
#
n_max = 1 
id_2 = np.eye(2) 
limit = 3 # maximum exponent -> 2^limit 
max_iterations = 100 # number of iterations
vec_gs = np.zeros(len(h_grid), dtype=float)
#
for ii,hh in enumerate(h_grid):
    (Ham_L,Ham_sing_L,Ham_R,Ham_sing_R, id_m, id_m_1,
                A_left, B_left,A_right, B_right) = initialize_operators(hh, sigma_x, sigma_z)
        # ______________________________________________
    for n_it in range(1,max_iterations):
        # before the enlargement 
        Ham_LR_int = np.kron(sigma_x,sigma_x)
        Ham_LR_int = np.kron(np.kron(A_left,Ham_LR_int),A_right)
        # enlargement 
        Ham_L = np.kron(Ham_L,id_2) + np.kron(A_left,Ham_sing_L) + np.kron(B_left,sigma_x)
        Ham_R = np.kron(id_2,Ham_R) + np.kron(Ham_sing_R,A_right) + np.kron(sigma_x,B_right)
        # ---- these will be used in the next iteration ----
        # 
        dim_half_sys = Ham_L.shape[0]
        id_half = np.eye(dim_half_sys)
        # update lef
        B_left = np.kron(A_left,sigma_x)
        A_left  = np.kron(A_left,id_2)
        # update right
        B_right = np.kron(sigma_x,A_right)
        A_right  = np.kron(id_2,A_right)
        # Total Hamiltonian 
        Ham_total = np.kron(Ham_L,id_half) + np.kron(id_half,Ham_R) + Ham_LR_int
        # Construction of the density matrix 
        eigenvals, eigenstates = LA.eigh(Ham_total)
        matrix_rho_AB = np.outer(eigenstates[:,0],eigenstates[:,0].conj())
        # reshape 
        new_shape = [dim_half_sys] * 4 
        tensor_rho_AB = matrix_rho_AB.reshape(new_shape)
        # trace out 
        matrix_rho_A = np.trace(tensor_rho_AB, axis1=1, axis2=3)
        # diagonalize the reduced density matrix 
        r_egvals, r_egvects = LA.eigh(matrix_rho_A)
        sorted_indices = np.argsort(r_egvals) 
        sorted_indices = sorted_indices[::-1]
        r_egvals = r_egvals[sorted_indices]
        r_egvects = r_egvects[:, sorted_indices]  
        # Choose how many df to keep 
        proj_states = 2 ** limit 
        Pr = r_egvects[:,:proj_states]
        Pr_T = np.transpose(Pr)
        Pr_H = np.conj(Pr_T)
        # projection 
        Ham_L = Pr_H @ ( Ham_L @ Pr)
        A_left = Pr_H @ ( A_left @ Pr)
        B_left = Pr_H @ ( B_left @ Pr)
        # 
        Ham_R = Pr_H @ ( Ham_R @ Pr)
        A_right = Pr_H @ ( A_right @ Pr)
        B_right = Pr_H @ ( B_right @ Pr)
    #
    hs_eigenvals,hs_eigenstates = LA.eigh(Ham_L)
    vec_gs[ii] = hs_eigenvals[0]
#
np.savetxt("GS_DMRG.txt",vec_gs)

