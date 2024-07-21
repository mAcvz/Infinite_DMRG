

import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as LA
import os

def generate_hamiltonian(N_max, h):
    """
    Generate the Hamiltonian for a quantum spin chain.

    Args:
    - N_max (int): Maximum number of spins in the chain.
    - h (float): Field strength.

    Returns:
    - Hamiltonian (2D array): The Hamiltonian matrix.
    """
    # Pauli matrices
    sigma_x = np.array([[0,1],[1,0]])
    sigma_z = np.array([[1,0],[0,-1]])
    identity = np.eye(2)
    
    list_H_o = [] 
    list_H_int = [] 

    for ii in range(1, N_max):
        # term on the left side - identities
        left_vec_id = np.eye(2 ** (ii-1))
        # single interaction term
        H_o = np.kron(left_vec_id, sigma_z)
        # add identities to the right
        for _ in range(ii, N_max):
            H_o = np.kron(H_o, identity)
        # store the field term
        list_H_o.append(H_o)
        # term on the left side - identities
        H_int = np.kron(left_vec_id, sigma_x)
        # insert the interaction term
        H_int = np.kron(H_int, sigma_x)
        # add identities to the right
        for _ in range(ii+1, N_max):
            H_int = np.kron(H_int, identity)
        # store the interaction term
        list_H_int.append(H_int)

    # Add the last term on the right side
    max_left_vec_id = np.eye(2 ** (N_max-1))
    H_o = np.kron(max_left_vec_id, sigma_z)
    list_H_o.append(H_o)
    # multiply by the constant field 

    # convert to numpy array
    vec_H_o = np.array(list_H_o)
    vec_H_o = vec_H_o * h
    vec_H_int = np.array(list_H_int)
    # Element wise sum 
    Hamiltonian = np.sum(vec_H_o, axis=0) + np.sum(vec_H_int, axis=0)
    
    return Hamiltonian
