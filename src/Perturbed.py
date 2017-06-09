from src import Steady
from src.util.Util import *
from src.util.MatrixPerUtil import *
import matplotlib.pyplot as plt

class Perturbed:

    def __init__(self, steadyInit = 0):
        print('Setting up Perturbed...', flush=True)
        self.steady = Steady.Steady(steadyInit)
        self.steady.run()
        print('\n...done')

    def run(self, loops = numLoops, verb='n'):
        print()
        if verb == 'v': 
            print("Running Perturbed...", flush=True)
        else:
            print("Running Perturbed...", end='', flush=True)

        for i in range(1, loops + 1):

            if verb == 'v': print("\tIteration number %i" % i)

            if verb == 'v': print("\t\tSolving for dY...") 
            del_dI = getdI()

            if verb == 'v': print("\t\tAdding dY to dI...") 
            dI[1:-1] += del_dI

    def getdI(self, indRange = 'i'):
        if indRange == 'i':
            return dI[1:-1]
        elif indRange == 'ip1':
            return dI[2:]
        elif indRange == 'im1':
            return dI[:-2]
        else:
            print("Index range error in I")
            exit(1)

    def plotdI_z(self, ir = 0, j = 0, k = 'all'):
        if k == "all":
            plt.plot(zGrid, dI[ir, 1:-1, j, :])
        elif j =="all":
            plt.plot(zGrid, dI[ir, 1:-1, :, k])
        else:
            plt.plot(zGrid, dI[ir, 1:-1, j, k])

    def render(self):
        print()
        plt.show(block=False)
        input('\tPress <ENTER> to continue')
        print()
