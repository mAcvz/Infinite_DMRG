
import numpy as np

def initialize_operators(hh, sigma_x, sigma_z):
    Ham_L = hh * sigma_z
    Ham_sing_L = hh * sigma_z
    Ham_R = hh * sigma_z
    Ham_sing_R = hh * sigma_z
    
    id_m = np.eye(2)
    id_m_1 = 1
    
    A_left = np.eye(2)
    B_left = sigma_x
    
    A_right = np.eye(2)
    B_right = sigma_x
    
    return Ham_L, Ham_sing_L, Ham_R, Ham_sing_R, id_m, id_m_1, A_left, B_left, A_right, B_right

