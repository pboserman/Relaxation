from src.util.MatrixStUtil import *
import matplotlib.pyplot as plt

from warnings import filterwarnings
filterwarnings("ignore")

# STEADY OBJECT CLASS TO RUN AND PLOT STEADY RESULTS

class Steady:

    def __init__(self, initial = 0):
        I[:] = initial

	# RUN PROGRAM FOR loops ITERATIONS
    def run(self, loops = numLoops, verb = 'n'):
        print()
        if verb == 'v': 
            print("Running Steady...", flush=True)
        else:
            print("Running Steady...", end='', flush=True)

        for i in range(1, loops + 1):

            if verb == 'v': print("\tIteration number %i" % i)

            if verb == 'v': print("\t\tSolving for dY...") 
            del_I = getdI()

            if verb == 'v': print("\t\tAdding dY to I...") 
            I[1:-1] += del_I

        if verb != 'v':
            print('done', flush = True)

    def getI(self, indRange = 'i'):
        if indRange == 'i':
            return I[1:-1]
        elif indRange == 'ip1':
            return I[2:]
        elif indRange == 'im1':
            return I[:-2]
        else:
            print("Index range error in I")
            exit(1)

    def plotI_z(self, j = 0, k = 'all'):
        if k == "all":
            plt.plot(zGrid, I[1:-1, j, :])
        elif j =="all":
            plt.plot(zGrid, I[1:-1, :, k])
        else:
            plt.plot(zGrid, I[1:-1, j, k])

    def render(self):
        print()
        plt.show(block=False)
        input('\tPress <ENTER> to continue')
        print()




