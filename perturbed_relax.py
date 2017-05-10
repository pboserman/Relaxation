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

I = np.zeros((zSize+2, xSize, muSize))

Xi = 0
# Theta_0 = 0
omega = 0.6
tau = 5
# flip = 1 # Turn Compton on and off

delZ = zGrid[1]-zGrid[0]
delX = xGrid[1]-xGrid[0]

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

def n_e(z):
    return (1-z**2)**3

def dn_e(z):
    return (1-z**2)**2 * (1-7*z**2)

def Temp(z):
    return (1-z**2)

def eta(z, x, which):
    if which == 'r':
        # delta eta term
        res = 1/(3)*(1-z**2)**3.5*(1-7*z**2)*(5.5*(1-z**2) + np.power(10, x))
    else:
        # eta term in imaginary part
        res = (1-z**2)**4.5*(2*(1-z**2) + np.power(10, x))
    return np.nan_to_num(res*np.exp(np.power(10, x)/(1-z**2)))

def weight(mu):
    muVals, weightVals = legpoly.legendre.leggauss(muSize)
    return weightVals

def G(z, mu):
    return -tau*n_e(z)/mu

# def A(z, x, mu, mu_p):
#     if len(z.shape) != 4:
#         z = z.reshape(-1, 1, 1, 1)
#     if len(x.shape) != 4:
#         x = x.reshape(1, -1, 1, 1)
#     if len(mu.shape) != 4:
#         mu = mu.reshape(1, 1, -1, 1)
#     if len(mu_p.shape) != 4:
#         mu_p = mu_p.reshape(1, 1, 1, -1)

#     # print(x.shape, z.shape, mu.shape)
#     return tau*n_e(z)/mu*(3./16.)*(alpha(mu, mu_p)*weight(mu_p))

def S_nu(mu, I_p):
    muGridp_4 = muGrid.reshape(1, 1, 1, -1)
    return 3/16*alpha(mu, muGridp_4)*I_p*weight(muGridp_4)

def F(z, x, mu, which):
    if which == 'r':
        return eta(z, x, which)/mu
    else:
        return omega*z*eta(z, x, which)

################################################

zGrid_3 = zGrid.reshape(-1, 1, 1)
xGrid_3 = xGrid.reshape(1, -1, 1)
muGrid_3 = muGrid.reshape(1, 1, -1)

zGrid_4 = zGrid.reshape(-1, 1, 1, 1)
xGrid_4 = xGrid.reshape(1, -1, 1, 1)
muGrid_4 = muGrid.reshape(1, 1, -1, 1)
muGridp_4 = muGrid.reshape(1, 1, 1, -1)







