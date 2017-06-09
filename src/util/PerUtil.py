import numpy as np
from opt import *
from src.util.Util import *
from numpy import newaxis as new


# ANGLE FUNCTIONS
########################################################################################################

def gamma(mu, mu_p):
    # Velocity I term
    return 8*mu-mu_p-5*mu**2*mu_p-8*mu_p**2*mu-4*mu**3-mu_p**3+12*mu**3*mu_p**2+3*mu_p**3*mu**2

def delta(mu, mu_p):
    # Velocity dI/dx term
    return 3*(mu_p-mu)+mu**3-mu_p**3+mu_p**2*mu-mu**2*mu_p+3*mu_p**3*mu**2-3*mu**3*mu_p**2

# PERTURBED SOURCE FUNCTIONS
########################################################################################################

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


