import numpy as np
import matplotlib.pyplot as plt

def plotBlaesFile(filename, axis = 'z'):
    fileArray = np.loadtxt(filename).T
    muVals = np.unique(fileArray[0])
    xVals = np.unique(fileArray[1])
    zVals = np.unique(fileArray[2])

    if axis == 'z':
        print(xVals)
        x = input("Enter x value from list (all for all x): ")
    if axis == 'x':
        print(zVals)
        z = input("Enter z value from list (all for all z): ")
    print(muVals)
    mu = input("Enter mu value from list ('all' for all mu): ")

    if axis == 'z':
        x = float(x) if x != 'all' else x
    else:
        z = float(z) if z != 'all' else z
    mu = float(mu) if mu != 'all' else mu

    if axis == 'x':
        fileArray[1] = np.log10(fileArray[1])
        # fileArray[3, fileArray[0] > 0] = np.log10(fileArray[3, fileArray[0] > 0])
        plotStyle = 'bo'
    else:
        plotStyle = 'b.'



    if axis == 'z' and x == 'all':
        plotIndexes = np.where(fileArray[0] == mu)
        plt.plot(fileArray[2, plotIndexes], fileArray[3, plotIndexes], plotStyle)

    elif axis == 'x' and z == 'all':
        plotIndexes = np.where(fileArray[2] == z)
        plt.plot(fileArray[1, plotIndexes], fileArray[3, plotIndexes], plotStyle)

    elif mu == 'all':
        if axis == 'z':
            plotIndexes = np.where(fileArray[1] == x)
            plt.plot(fileArray[2, plotIndexes], fileArray[3, plotIndexes], plotStyle)
        else:
            plotIndexes = np.where( (fileArray[2] == z) )
            plt.plot(fileArray[1, plotIndexes], fileArray[3, plotIndexes], plotStyle)

    else:
        if axis == 'z':
            plotIndexes = np.where( (fileArray[1] == x) and (fileArray[0] == mu) )
            plt.plot(fileArray[2, plotIndexes], fileArray[3, plotIndexes], plotStyle)
        else:
            plotIndexes = np.where( (fileArray[2] == z) and (fileArray[0] == mu) )
            plt.plot(fileArray[1, plotIndexes], fileArray[3, plotIndexes], plotStyle)


