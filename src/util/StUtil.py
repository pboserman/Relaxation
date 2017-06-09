import numpy as np
from opt import *
from src.util.Util import *

# BASIC UTILITY FUNCTIONS FOR STEADY

def Temp(z):
    return (1-z**2)

def eta(z, x):
    # T = Temp(z)
    # res = n_e(z)**2/np.sqrt(T)*np.exp(-np.power(10, x)/T)
    res = (1-z**2)**(5.5)*np.exp(-np.power(10, x)/(1-z**2))
    return np.nan_to_num(res)

def G(z, mu):
    return tau*n_e(z)/mu

def A(z, x, mu, mu_p):
    if len(z.shape) != 4:
        z = z.reshape(-1, 1, 1, 1)
    if len(x.shape) != 4:
        x = x.reshape(1, -1, 1, 1)
    if len(mu.shape) != 4:
        mu = mu.reshape(1, 1, -1, 1)
    if len(mu_p.shape) != 4:
        mu_p = mu_p.reshape(1, 1, 1, -1)

    return tau*n_e(z)/mu*(3./16.)*(alpha(mu, mu_p)*weight(mu_p))


def F(z, x, mu):
    return eta(z, x)/mu