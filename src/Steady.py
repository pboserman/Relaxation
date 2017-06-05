
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

zGrid_3 = zGrid.reshape(-1, 1, 1)
xGrid_3 = xGrid.reshape(1, -1, 1)
muGrid_3 = muGrid.reshape(1, 1, -1)

zGrid_4 = zGrid.reshape(-1, 1, 1, 1)
xGrid_4 = xGrid.reshape(1, -1, 1, 1)
muGrid_4 = muGrid.reshape(1, 1, -1, 1)
muGridp_4 = muGrid.reshape(1, 1, 1, -1)

def getI(i, j, k):
    if i > zSize - 1:
        if muGrid[k] < 0:
            return 0
        else:
            print("Shouldn't be called")
            exit(1)
    if i < 0:
        if muGrid[k] > 0:
            return 0
        else:
            print("Shouldn't be called")
            exit(1)
    return I[i][j][k]

def getI(indRange = 'i'):
    if indRange == 'i':
        return I[1:-1]
    elif indRange == 'ip1':
        return I[2:]
    elif indRange == 'im1':
        return I[:-2]
    else:
        print("Index range error in I")
        exit(1)

def plotI_z(j, k):
    if k == "all":
        plt.plot(zGrid, I[1:-1, j, :])
    elif j =="all":
        plt.plot(zGrid, I[1:-1, :, k])
    else:
        plt.plot(zGrid, I[1:-1, j, k])

Id_x = np.eye(xSize, xSize).reshape(1, xSize, xSize, 1, 1)
Id_zp1 = np.eye(zSize, zSize, 1).reshape(zSize, zSize, 1, 1, 1, 1)
Id_zm1 = np.eye(zSize, zSize, -1).reshape(zSize, zSize, 1, 1, 1, 1)
Id_z = np.eye(zSize, zSize).reshape(zSize, zSize, 1, 1, 1, 1)

# ITERATION FOR SOLUTION
def run(loops = numLoops):
    for i in range(1, loops+1):

        print("Iteration number %i" % i)

        # S Matrix creation

        print("\tCreating S_M...")

        S_M = ( ( S_km1(zGrid_4, xGrid_4, muGrid_4, muGridp_4) * (muGrid_4 > 0) )[:, :, new] * Id_x)[:, new] * Id_zm1 # * (muGrid_4 > 0) 

        S_M += ( ( S_k(zGrid_4, xGrid_4, muGrid_4, muGridp_4, 'n') * (muGrid_4 > 0) )[:, :, new] * Id_x)[:, new] * Id_z # * (muGrid_4 > 0)

        S_M += ( ( S_k(zGrid_4, xGrid_4, muGrid_4, muGridp_4, 'p') * (muGrid_4 < 0) )[:, :, new] * Id_x)[:, new] * Id_z

        S_M += ( ( S_kp1(zGrid_4, xGrid_4, muGrid_4, muGridp_4) * (muGrid_4 < 0) )[:, :, new] * Id_x)[:, new] * Id_zp1

        S_M = S_M.swapaxes(1, 2).swapaxes(2, 4).swapaxes(3, 4).reshape(zSize*xSize*muSize, -1)


        print("\tSolving for dY...")

        dI = np.linalg.solve(S_M, -E_k(zGrid_3, xGrid_3, muGrid_3).reshape(-1)).reshape(zSize, xSize, muSize)


        print("\tAdding dY to I...")

        I[1:-1] += dI


##################################################################################################################