from src.util.StUtil import *
from src.util.Util import *
from opt import *
import numpy as np
from numpy import newaxis as new
# MATRIX UTILITY FUNCTIONS FOR STEADY

def E_k(z, x, mu):
    I_i = I[1:-1]
    I_ip1 = I[2:]
    I_im1 = I[:-2]
    
    E_1 = I_ip1 - I_i + delZ*G(z+0.5*delZ, mu).repeat(xSize, axis = 1)/2*(I_ip1+I_i) - delZ*F(z+0.5*delZ, x, mu) - delZ*np.sum(A(z, x, mu, muGridp_4)/2*np.expand_dims(I_i+I_ip1, axis = -2), axis = -1)
    E_2 = I_i - I_im1 + delZ*G(z-0.5*delZ, mu).repeat(xSize, axis = 1)/2*(I_im1+I_i) - delZ*F(z-0.5*delZ, x, mu) - delZ*np.sum(A(z, x, mu, muGridp_4)/2*np.expand_dims(I_i+I_im1, axis = -2), axis = -1)
    
    return E_1*(muGrid < 0) + E_2*(muGrid > 0)

def S_k(z, x, mu, mu_p, which):
    if which == 'p':
        return (G(z+0.5*delZ, mu).repeat(xSize, axis = 1)/2*delZ - 1)*(mu==mu_p) - A(z+0.5*delZ, x, mu, mu_p)/2*delZ
    if which == 'n':
        return (G(z-0.5*delZ, mu).repeat(xSize, axis = 1)/2*delZ + 1)*(mu==mu_p) - A(z-0.5*delZ, x, mu, mu_p)/2*delZ
    else:
        print("Option not available")
        exit(1)

def S_km1(z, x, mu, mu_p):
    return (G(z-0.5*delZ, mu).repeat(xSize, axis = 1)/2*delZ - 1)*(mu==mu_p) - A(z-0.5*delZ, x, mu, mu_p)/2*delZ

def S_kp1(z, x, mu, mu_p):
    return (G(z+0.5*delZ, mu).repeat(xSize, axis = 1)/2*delZ + 1)*(mu==mu_p) - A(z+0.5*delZ, x, mu, mu_p)/2*delZ

def make_S_M():
    S_M = ( ( S_km1(zGrid_4, xGrid_4, muGrid_4, muGridp_4) * (muGrid_4 > 0) )[:, :, new] * Ids_x)[:, new] * Ids_zm1
    
    S_M += ( ( S_k(zGrid_4, xGrid_4, muGrid_4, muGridp_4, 'n') * (muGrid_4 > 0) )[:, :, new] * Ids_x)[:, new] * Ids_z 
    
    S_M += ( ( S_k(zGrid_4, xGrid_4, muGrid_4, muGridp_4, 'p') * (muGrid_4 < 0) )[:, :, new] * Ids_x)[:, new] * Ids_z
    
    S_M += ( ( S_kp1(zGrid_4, xGrid_4, muGrid_4, muGridp_4) * (muGrid_4 < 0) )[:, :, new] * Ids_x)[:, new] * Ids_zp1
    
    return S_M.swapaxes(1, 2).swapaxes(2, 4).swapaxes(3, 4).reshape(zSize*xSize*muSize, -1)

def getdI():
    S_M = make_S_M()
    del_I = np.linalg.solve(S_M, -E_k(zGrid_3, xGrid_3, muGrid_3).reshape(-1)).reshape(zSize, xSize, muSize)
    return del_I

