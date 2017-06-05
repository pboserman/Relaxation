import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
from scipy.integrate import quadrature as integrate
from scipy.integrate import ode
import time
import numpy.polynomial as legpoly
from warnings import filterwarnings
from numpy import newaxis as new

import relaxation as sted
sted.run(1)
sted.I[sted.I < 0] = 0

# filterwarnings("ignore")


# x = hv/kT

def getMuGrid(mu_size):
    muVals, weightVals = legpoly.legendre.leggauss(mu_size)
    return muVals

numLoops = 1

zSize = 100
xSize = 10
muSize = 8

x0 = -1
xMax = 1

z0 = -1
zMax = 1

zGrid = np.linspace(z0, zMax, zSize)
xGrid = np.linspace(x0, xMax, xSize)
muGrid = getMuGrid(muSize)

I = sted.I
dI = np.zeros((2, zSize+2, xSize, muSize))

Xi = 0
# Theta_0 = 0
omega = 0.6
tau = 5
# flip = 1 # Turn Compton on and off

delZ = zGrid[1]-zGrid[0]
delX = xGrid[1]-xGrid[0]

# ANGLE FUNCTIONS
########################################################################################################
def alpha(mu, mu_p):
    return 3*(mu*mu_p)**2-mu**2-mu_p**2+3

def beta(mu, mu_p):
    return 5*(mu*mu_p+(mu*mu_p)**3)-3*(mu**3*mu_p-mu*mu_p**3)

def epsilon(mu, mu_p):
    return -3*alpha(mu, mu_p)+8+2*beta(mu, mu_p)-8*mu*mu_p

def gamma(mu, mu_p):
    # Velocity I term
    return 8*mu-mu_p-5*mu**2*mu_p-8*mu_p**2*mu-4*mu**3-mu_p**3+12*mu**3*mu_p**2+3*mu_p**3*mu**2

def delta(mu, mu_p):
    # Velocity dI/dx term
    return 3*(mu_p-mu)+mu**3-mu_p**3+mu_p**2*mu-mu**2*mu_p+3*mu_p**3*mu**2-3*mu**3*mu_p**2

# EMISSIVITY, DENSITY, VELOCITY, etc... FUNCTIONS
########################################################################################################

def n_e(z):
    return (1-z**2)**3

def dn_e(z):
    return (1-z**2)**2 * (1-7*z**2)

def eta(z, x, which):
    if which == 0:
        # delta eta term
        res = 1/(3)*(1-z**2)**3.5*(1-7*z**2)*(5.5*(1-z**2) + np.power(10, x))
    else:
        # eta term in imaginary part
        res = (1-z**2)**4.5*(2*(1-z**2) + np.power(10, x))
    return np.nan_to_num(res*np.exp(np.power(10, x)/(1-z**2)))

def dV(z):
    return omega*z

def weight(mu):
    muVals, weightVals = legpoly.legendre.leggauss(muSize)
    return weightVals

# OTHER HELPER/WRAPPER FUNCTIONS 
########################################################################################################

def getI_k(I_arr, offset = 0):
    if offset == 0:
        return I_arr
    if offset == 1:
        return np.append(I_arr, (2*I_arr[:, -1] - I_arr[:, -2]).reshape(zSize, 1, muSize), axis = 1)[:, 1:]
    if offset == -1:
        return np.append((2*I_arr[:, 0] - I_arr[:, 1]).reshape(zSize, 1, muSize), I_arr, axis = 1)[:, :-1]


def S_nu(I_p):
    muGrid_4 = muGrid.reshape(1, 1, -1, 1)
    muGridp_4 = muGrid.reshape(1, 1, 1, -1)
    res = 3/16*alpha(muGrid_4, muGridp_4)*I_p[:, :, new, :]*weight(muGridp_4)
    return np.sum(res, axis = -1)

def dS_nu(dI_p, which = 0, I_p2 = [], I_p3 = []):
    zGrid_4 = zGrid.reshape(-1, 1, 1, 1)
    muGrid_4 = muGrid.reshape(1, 1, -1, 1)
    muGridp_4 = muGrid.reshape(1, 1, 1, -1)
    res = 3/16*alpha(muGrid_4, muGridp_4)*dI_p[:, :, new, :]*weight(muGridp_4)
    if which == 1:
        res += 3/16*dV(zGrid_4)*weight(muGridp_4)*(gamma(muGrid_4, muGridp_4)*I_p2[:, :, new, :] + delta(muGrid_4, muGridp_4)*I_p3[:, :, new, :])
    return np.sum(res, axis = -1)

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

zGrid_3 = zGrid.reshape(-1, 1, 1)
xGrid_3 = xGrid.reshape(1, -1, 1)
muGrid_3 = muGrid.reshape(1, 1, -1)

zGrid_4 = zGrid.reshape(-1, 1, 1, 1)
xGrid_4 = xGrid.reshape(1, -1, 1, 1)
muGrid_4 = muGrid.reshape(1, 1, -1, 1)
muGridp_4 = muGrid.reshape(1, 1, 1, -1)

# ID Creations

Id_x = np.eye(xSize, xSize).reshape(1, 1, xSize, xSize, 1, 1)
Id_zp1 = np.eye(zSize, zSize, 1).reshape(1, zSize, zSize, 1, 1, 1, 1)
Id_zm1 = np.eye(zSize, zSize, -1).reshape(1, zSize, zSize, 1, 1, 1, 1)
Id_z = np.eye(zSize, zSize).reshape(1, zSize, zSize, 1, 1, 1, 1)
Id_ir = np.eye(2, 2).reshape(2, 2, 1, 1, 1, 1, 1, 1)
Idd_ir = (np.eye(2, 2, 1) + np.eye(2, 2, -1)).reshape(2, 2, 1, 1, 1, 1, 1, 1)

def run(loops = numLoops):
    for i in range(1, loops+1):

        print("Iteration number %i" % i)

        # S Matrix creation

        print("\tCreating S_M...")

        S_M = ( ( ( S_km1_diag(zGrid_4, xGrid_4, muGrid_4, muGridp_4) * (muGrid_4 > 0) )[:, :, :, new] * Id_x )[:, :, new] * Id_zm1 )[:, new] * Id_ir

        S_M += ( ( ( S_k_diag(zGrid_4, xGrid_4, muGrid_4, muGridp_4, 'n') * (muGrid_4 > 0) )[:, :, :, new] * Id_x )[:, :, new] * Id_z )[:, new] * Id_ir

        S_M += ( ( ( S_k_diag(zGrid_4, xGrid_4, muGrid_4, muGridp_4, 'p') * (muGrid_4 < 0) )[:, :, :, new] * Id_x )[:, :, new] * Id_z )[:, new] * Id_ir

        S_M += ( ( ( S_kp1_diag(zGrid_4, xGrid_4, muGrid_4, muGridp_4) * (muGrid_4 < 0) )[:, :, :, new] * Id_x )[:, :, new] * Id_zp1 )[:, new] * Id_ir


        S_M += ( ( ( S_km1_nondiag(zGrid_4, xGrid_4, muGrid_4, muGridp_4) * (muGrid_4 > 0) )[:, :, :, new] * Id_x )[:, :, new] * Id_zm1 )[:, new] * Idd_ir

        S_M += ( ( ( S_k_nondiag(zGrid_4, xGrid_4, muGrid_4, muGridp_4, 'n') * (muGrid_4 > 0) )[:, :, :, new] * Id_x )[:, :, new] * Id_z )[:, new] * Idd_ir

        S_M += ( ( ( S_k_nondiag(zGrid_4, xGrid_4, muGrid_4, muGridp_4, 'p') * (muGrid_4 < 0) )[:, :, :, new] * Id_x )[:, :, new] * Id_z )[:, new] * Idd_ir

        S_M += ( ( ( S_kp1_nondiag(zGrid_4, xGrid_4, muGrid_4, muGridp_4) * (muGrid_4 < 0) )[:, :, :, new] * Id_x )[:, :, new] * Id_zp1 )[:, new] * Idd_ir


        # TODO: EDIT FOR CORRECT SWAPPING WITH 6 DIMENSIONS
        S_M = S_M.swapaxes(1, 2).swapaxes(2, 4).swapaxes(3, 6).swapaxes(6, 5).reshape(2*zSize*xSize*muSize, -1)

        print("\tSolving for dY..")

        ddI = np.linalg.solve(S_M, -E_k(zGrid_3, xGrid_3, muGrid_3).reshape(-1)).reshape(2, zSize, xSize, muSize)

        print("\tAdding dY to I...")
        
        dI[:, 1:-1] += ddI

def getdI(indRange = 'i'):
    if indRange == 'i':
        return dI[1:-1]
    elif indRange == 'ip1':
        return dI[2:]
    elif indRange == 'im1':
        return dI[:-2]
    else:
        print("Index range error in I")
        exit(1)

def plotdI_z(ir, j, k):
    if k == "all":
        plt.plot(zGrid, dI[ir, 1:-1, j, :])
    elif j =="all":
        plt.plot(zGrid, dI[ir, 1:-1, :, k])
    else:
        plt.plot(zGrid, dI[ir, 1:-1, j, k])
