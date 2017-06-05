import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
from scipy.integrate import quadrature as integrate
from scipy.integrate import ode
import time

import numpy.polynomial as legpoly

from warnings import filterwarnings

from numpy import newaxis as new

filterwarnings("ignore")


# x = hv/kT

# def getMuGrid(mu_size):
#     muVals, weightVals = legpoly.legendre.leggauss(mu_size)
#     return muVals

# numLoops = 1

# zSize = 100
# xSize = 10
# muSize = 8

# x0 = -1
# xMax = 1

# z0 = -1
# zMax = 1

# zGrid = np.linspace(z0, zMax, zSize)
# xGrid = np.linspace(x0, xMax, xSize)
# muGrid = getMuGrid(muSize)

# I = np.zeros((zSize+2, xSize, muSize))

# Xi = 0
# # Theta_0 = 0
# omega = 0.6
# tau = 5
# # flip = 1 # Turn Compton on and off

# delZ = zGrid[1]-zGrid[0]
# delX = xGrid[1]-xGrid[0]

# def getI(i, j, k):
#     if i > zSize - 1:
#         if muGrid[k] < 0:
#             return 0
#         else:
#             print("Shouldn't be called")
#             exit(1)
#     if i < 0:
#         if muGrid[k] > 0:
#             return 0
#         else:
#             print("Shouldn't be called")
#             exit(1)
#     return I[i][j][k]



# def alpha(mu, mu_p):
#     return 3*(mu*mu_p)**2-mu**2-mu_p**2+3

# def beta(mu, mu_p):
#     return 5*(mu*mu_p+(mu*mu_p)**3)-3*(mu**3*mu_p-mu*mu_p**3)

# def epsilon(mu, mu_p):
#     return -3*alpha(mu, mu_p)+8+2*beta(mu, mu_p)-8*mu*mu_p

# def n_e(z):
#     return (1-z**2)**3

# def Temp(z):
#     return (1-z**2)

# def eta(z, x):
#     # T = Temp(z)
#     # res = n_e(z)**2/np.sqrt(T)*np.exp(-np.power(10, x)/T)
#     res = (1-z**2)**(5.5)*np.exp(-np.power(10, x)/(1-z**2))
#     return np.nan_to_num(res)

# def weight(mu):
#     muVals, weightVals = legpoly.legendre.leggauss(muSize)
#     return weightVals

# def G(z, mu):
#     return tau*n_e(z)/mu

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



# def F(z, x, mu):
#     return eta(z, x)/mu

################################################

# zGrid_3 = zGrid.reshape(-1, 1, 1)
# xGrid_3 = xGrid.reshape(1, -1, 1)
# muGrid_3 = muGrid.reshape(1, 1, -1)

# zGrid_4 = zGrid.reshape(-1, 1, 1, 1)
# xGrid_4 = xGrid.reshape(1, -1, 1, 1)
# muGrid_4 = muGrid.reshape(1, 1, -1, 1)
# muGridp_4 = muGrid.reshape(1, 1, 1, -1)

# def E_k(z, x, mu):
#     I_i = I[1:-1]
#     I_ip1 = I[2:]
#     I_im1 = I[:-2]
    
#     E_1 = I_ip1 - I_i + delZ*G(z+0.5*delZ, mu).repeat(xSize, axis = 1)/2*(I_ip1+I_i) - delZ*F(z+0.5*delZ, x, mu) - delZ*np.sum(A(z, x, mu, muGridp_4)/2*np.expand_dims(I_i+I_ip1, axis = -2), axis = -1)
#     E_2 = I_i - I_im1 + delZ*G(z-0.5*delZ, mu).repeat(xSize, axis = 1)/2*(I_im1+I_i) - delZ*F(z-0.5*delZ, x, mu) - delZ*np.sum(A(z, x, mu, muGridp_4)/2*np.expand_dims(I_i+I_im1, axis = -2), axis = -1)
    
#     return E_1*(muGrid < 0) + E_2*(muGrid > 0)

# def S_k(z, x, mu, mu_p, which):
#     if which == 'p':
#         return (G(z+0.5*delZ, mu).repeat(xSize, axis = 1)/2*delZ - 1)*(mu==mu_p) - A(z+0.5*delZ, x, mu, mu_p)/2*delZ
#     if which == 'n':
#         return (G(z-0.5*delZ, mu).repeat(xSize, axis = 1)/2*delZ + 1)*(mu==mu_p) - A(z-0.5*delZ, x, mu, mu_p)/2*delZ
#     else:
#         print("Option not available")
#         exit(1)

# def S_km1(z, x, mu, mu_p):
#     return (G(z-0.5*delZ, mu).repeat(xSize, axis = 1)/2*delZ - 1)*(mu==mu_p) - A(z-0.5*delZ, x, mu, mu_p)/2*delZ

# def S_kp1(z, x, mu, mu_p):
#     return (G(z+0.5*delZ, mu).repeat(xSize, axis = 1)/2*delZ + 1)*(mu==mu_p) - A(z+0.5*delZ, x, mu, mu_p)/2*delZ


# ID Creations

# Id_x = np.eye(xSize, xSize).reshape(1, xSize, xSize, 1, 1)
# Id_zp1 = np.eye(zSize, zSize, 1).reshape(zSize, zSize, 1, 1, 1, 1)
# Id_zm1 = np.eye(zSize, zSize, -1).reshape(zSize, zSize, 1, 1, 1, 1)
# Id_z = np.eye(zSize, zSize).reshape(zSize, zSize, 1, 1, 1, 1)

# # ITERATION FOR SOLUTION
# def run(loops = numLoops):
#     for i in range(1, loops+1):

#         print("Iteration number %i" % i)

#         # S Matrix creation

#         print("\tCreating S_M...")

#         S_M = ( ( S_km1(zGrid_4, xGrid_4, muGrid_4, muGridp_4) * (muGrid_4 > 0) )[:, :, new] * Id_x)[:, new] * Id_zm1 # * (muGrid_4 > 0) 

#         S_M += ( ( S_k(zGrid_4, xGrid_4, muGrid_4, muGridp_4, 'n') * (muGrid_4 > 0) )[:, :, new] * Id_x)[:, new] * Id_z # * (muGrid_4 > 0)

#         S_M += ( ( S_k(zGrid_4, xGrid_4, muGrid_4, muGridp_4, 'p') * (muGrid_4 < 0) )[:, :, new] * Id_x)[:, new] * Id_z

#         S_M += ( ( S_kp1(zGrid_4, xGrid_4, muGrid_4, muGridp_4) * (muGrid_4 < 0) )[:, :, new] * Id_x)[:, new] * Id_zp1

#         S_M = S_M.swapaxes(1, 2).swapaxes(2, 4).swapaxes(3, 4).reshape(zSize*xSize*muSize, -1)


#         print("\tSolving for dY...")

#         dI = np.linalg.solve(S_M, -E_k(zGrid_3, xGrid_3, muGrid_3).reshape(-1)).reshape(zSize, xSize, muSize)


#         print("\tAdding dY to I...")

#         I[1:-1] += dI

# dYMax = np.amax(dY)

# maxI = I[np.where(dY == dYMax)]

# print(dY)

# print("Max dY is: %.2f" % ( dYMax ) )

# print(np.where(dY == dYMax))

# for j in range(xSize):
#     plt.figure("X = %.2f" % (10**xGrid[j]))
#     plt.plot(zGrid, I[1:-1, j, :])

# def getI(indRange = 'i'):
#     if indRange == 'i':
#         return I[1:-1]
#     elif indRange == 'ip1':
#         return I[2:]
#     elif indRange == 'im1':
#         return I[:-2]
#     else:
#         print("Index range error in I")
#         exit(1)

# def plotI_z(j, k):
#     if k == "all":
#         plt.plot(zGrid, I[1:-1, j, :])
#     elif j =="all":
#         plt.plot(zGrid, I[1:-1, :, k])
#     else:
#         plt.plot(zGrid, I[1:-1, j, k])





