import numpy.polynomial as legpoly
import numpy as np
from opt import *

# BASIC SHARED UTILITY FUNCTIONS FOR STEADY AND PERTURBED

def getMuGrid(mu_size):
    muVals, weightVals = legpoly.legendre.leggauss(mu_size)
    return muVals

def weight(mu):
    muVals, weightVals = legpoly.legendre.leggauss(muSize)
    return weightVals

def n_e(z):
    return (1-z**2)**3

def alpha(mu, mu_p):
    return 3*(mu*mu_p)**2-mu**2-mu_p**2+3

def beta(mu, mu_p):
    return 5*(mu*mu_p+(mu*mu_p)**3)-3*(mu**3*mu_p-mu*mu_p**3)

def epsilon(mu, mu_p):
    return -3*alpha(mu, mu_p)+8+2*beta(mu, mu_p)-8*mu*mu_p

# GRID CREATION
zGrid = np.linspace(z0, zMax, zSize)
xGrid = np.linspace(x0, xMax, xSize)
muGrid = getMuGrid(muSize)
I = np.zeros((zSize+2, xSize, muSize))
dI = np.zeros((2, zSize+2, xSize, muSize))

# DELTA VARS
delZ = zGrid[1]-zGrid[0]
delX = xGrid[1]-xGrid[0]



# OTHER UTILITY GRIDS AND MATRICES
Ids_x = np.eye(xSize, xSize).reshape(1, xSize, xSize, 1, 1)
Ids_zp1 = np.eye(zSize, zSize, 1).reshape(zSize, zSize, 1, 1, 1, 1)
Ids_zm1 = np.eye(zSize, zSize, -1).reshape(zSize, zSize, 1, 1, 1, 1)
Ids_z = np.eye(zSize, zSize).reshape(zSize, zSize, 1, 1, 1, 1)

Idp_x = np.eye(xSize, xSize).reshape(1, 1, xSize, xSize, 1, 1)
Idp_zp1 = np.eye(zSize, zSize, 1).reshape(1, zSize, zSize, 1, 1, 1, 1)
Idp_zm1 = np.eye(zSize, zSize, -1).reshape(1, zSize, zSize, 1, 1, 1, 1)
Idp_z = np.eye(zSize, zSize).reshape(1, zSize, zSize, 1, 1, 1, 1)
Idp_ir = np.eye(2, 2).reshape(2, 2, 1, 1, 1, 1, 1, 1)
Idp_ir_inv = (np.eye(2, 2, 1) + np.eye(2, 2, -1)).reshape(2, 2, 1, 1, 1, 1, 1, 1)

zGrid_3 = zGrid.reshape(-1, 1, 1)
xGrid_3 = xGrid.reshape(1, -1, 1)
muGrid_3 = muGrid.reshape(1, 1, -1)

zGrid_4 = zGrid.reshape(-1, 1, 1, 1)
xGrid_4 = xGrid.reshape(1, -1, 1, 1)
muGrid_4 = muGrid.reshape(1, 1, -1, 1)
muGridp_4 = muGrid.reshape(1, 1, 1, -1)