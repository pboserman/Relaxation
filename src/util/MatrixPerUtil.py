from src.util.PerUtil import *
from src.util.Util import *
from opt import *
import numpy as np
from numpy import newaxis as new

# RELAXATION METHOD FUNCTIONS
##########################################################################################################

def E_k(z, x, mu):
    dI_i = dI[:, 1:-1]
    dI_ip1 = dI[:, 2:]
    dI_im1 = dI[:, :-2]

    dI_avg_1 = (dI_i+dI_ip1)/2
    dI_avg_2 = (dI_i+dI_im1)/2

    I_i = I[1:-1]
    I_ip1 = I[2:]
    I_im1 = I[:-2]

    I_avg_1 = (I_i+I_ip1)/2
    I_avg_2 = (I_i+I_im1)/2

    dIdx_1 = getI_k(I_avg_1, 1)/2 - getI_k(I_avg_1, -1)/2
    dIdx_2 = getI_k(I_avg_2, 1)/2 - getI_k(I_avg_2, -1)/2

    E_1_r = -omega*(dI_avg_1[1]) + tau*dn_e(z)*(-I_avg_1 + S_nu(I_avg_1)) + tau*n_e(z)*(-dI_avg_1[0] + dS_nu(dI_avg_1[0], 0)) + eta(z, x, 0)

    E_2_r = -omega*(dI_avg_2[1]) + tau*dn_e(z)*(-I_avg_2 + S_nu(I_avg_2)) + tau*n_e(z)*(-dI_avg_2[0] + dS_nu(dI_avg_2[0], 0)) + eta(z, x, 0)

    E_1_i = omega*(dI_avg_1[0]) + tau*n_e(z)*(-dI_avg_1[1] + dS_nu(dI_avg_1[1], 1, I_avg_1, x/delX*dIdx_1) + mu*dV(z)*I_avg_1) + mu*dV(z)*eta(z, x, 1)

    E_2_i = omega*(dI_avg_2[0]) + tau*n_e(z)*(-dI_avg_2[1] + dS_nu(dI_avg_2[1], 1, I_avg_2, x/delX*dIdx_2) + mu*dV(z)*I_avg_2) + mu*dV(z)*eta(z, x, 1)

    E_1_r = dI_ip1[0] - dI_i[0] - delZ/mu*E_1_r
    E_2_r = dI_i[0] - dI_im1[0] - delZ/mu*E_2_r
    E_1_i = dI_ip1[1] - dI_i[1] - delZ/mu*E_1_i
    E_2_i = dI_i[1] - dI_im1[1] - delZ/mu*E_2_i

    return np.append(E_1_r[new], E_1_i[new], axis = 0) * (muGrid < 0) + np.append(E_2_r[new], E_2_i[new], axis = 0) * (muGrid > 0)

def S_k_diag(z, x, mu, mu_p, which):
    # Partial derivative with respect to dI_i

    res_r = (-tau*n_e(z)) * (mu == mu_p) + tau*n_e(z)*3/16*alpha(mu, mu_p)*weight(mu_p)

    if which == 'p':
        res_r = - 1 - delZ/mu*res_r/2
    else:
        res_r = 1 - delZ/mu*res_r/2

    res_i = res_r

    return np.append(res_r[new], res_i[new], axis = 0)

def S_km1_diag(z, x, mu, mu_p):
    # Partial derivative with respect to dI_ip1
    res_r = (-tau*n_e(z)) * (mu == mu_p) + tau*n_e(z)*3/16*alpha(mu, mu_p)*weight(mu_p)

    res_r = - 1 - delZ/mu*res_r/2

    res_i = res_r

    return np.append(res_r[new], res_i[new], axis = 0)

def S_kp1_diag(z, x, mu, mu_p):
    # Partial derivative with respect to other dI_im1
    res_r = (-tau*n_e(z)) * (mu == mu_p) + tau*n_e(z)*3/16*alpha(mu, mu_p)*weight(mu_p)

    res_r = 1 - delZ/mu*res_r/2

    res_i = res_r

    return np.append(res_r[new], res_i[new], axis = 0)

def S_k_nondiag(z, x, mu, mu_p, which):
    # Partial derivative with respect to other dI_i
    res_r = omega/2*delZ/mu
    res_i = -omega/2*delZ/mu

    return np.append(res_r[new], res_i[new], axis = 0)

def S_km1_nondiag(z, x, mu, mu_p):
    # Partial derivative with respect to other dI_ip1
    res_r = omega/2*delZ/mu
    res_i = -omega/2*delZ/mu

    return np.append(res_r[new], res_i[new], axis = 0)

def S_kp1_nondiag(z, x, mu, mu_p):
    # Partial derivative with respect to other dI_im1
    res_r = omega/2*delZ/mu
    res_i = -omega/2*delZ/mu

    return np.append(res_r[new], res_i[new], axis = 0)

################################################################################################################################################

def make_S_M():

    S_M = ( ( ( S_km1_diag(zGrid_4, xGrid_4, muGrid_4, muGridp_4) * (muGrid_4 > 0) )[:, :, :, new] * Idp_x )[:, :, new] * Idp_zm1 )[:, new] * Idp_ir
    S_M += ( ( ( S_k_diag(zGrid_4, xGrid_4, muGrid_4, muGridp_4, 'n') * (muGrid_4 > 0) )[:, :, :, new] * Idp_x )[:, :, new] * Idp_z )[:, new] * Idp_ir
    S_M += ( ( ( S_k_diag(zGrid_4, xGrid_4, muGrid_4, muGridp_4, 'p') * (muGrid_4 < 0) )[:, :, :, new] * Idp_x )[:, :, new] * Idp_z )[:, new] * Idp_ir
    S_M += ( ( ( S_kp1_diag(zGrid_4, xGrid_4, muGrid_4, muGridp_4) * (muGrid_4 < 0) )[:, :, :, new] * Idp_x )[:, :, new] * Idp_zp1 )[:, new] * Idp_ir
    S_M += ( ( ( S_km1_nondiag(zGrid_4, xGrid_4, muGrid_4, muGridp_4) * (muGrid_4 > 0) )[:, :, :, new] * Idp_x )[:, :, new] * Idp_zm1 )[:, new] * Idp_ir_inv
    S_M += ( ( ( S_k_nondiag(zGrid_4, xGrid_4, muGrid_4, muGridp_4, 'n') * (muGrid_4 > 0) )[:, :, :, new] * Idp_x )[:, :, new] * Idp_z )[:, new] * Idp_ir_inv
    S_M += ( ( ( S_k_nondiag(zGrid_4, xGrid_4, muGrid_4, muGridp_4, 'p') * (muGrid_4 < 0) )[:, :, :, new] * Idp_x )[:, :, new] * Idp_z )[:, new] * Idp_ir_inv
    S_M += ( ( ( S_kp1_nondiag(zGrid_4, xGrid_4, muGrid_4, muGridp_4) * (muGrid_4 < 0) )[:, :, :, new] * Idp_x )[:, :, new] * Idp_zp1 )[:, new] * Idp_ir_inv

    S_M = S_M.swapaxes(1, 2).swapaxes(2, 4).swapaxes(3, 6).swapaxes(6, 5).reshape(2*zSize*xSize*muSize, -1)
    return S_M

def getdI():
    S_M = make_S_M()
    del_dI = np.linalg.solve(S_M, -E_k(zGrid_3, xGrid_3, muGrid_3).reshape(-1)).reshape(2, zSize, xSize, muSize)
    return del_dI
