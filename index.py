# Handles App object creation and Plotting
from src import Steady
from src import Perturbed

""" *****************************************
USAGE
-----

* IN ORDER TO CREATE A NEW INSTANCE OF STEADY/PERTURBED

    steady = Steady.Steady(initValue)
    perturbed = Perturbed.Perturbed(initValue)


    initValue -> initial array/number guess 
                 to be set for steady solution
    
    ** if not specified it defaults to 0 **


* IN ORDER TO RUN

    steady.run(loops = 1, verb = 'n')
    perturbed.run(loops = 1, verb = 'n')


    loops -> Number of iterations for relaxation method 
    verb -> Verbose mode: 'v' for verbose

    ** loops defaults to numLoops from opt.py **
    ** verb defaults to non-verbose **

* OTHER AVAILABLE METHODS

    steady.plotI_z(j, k)
    perturbed.plotdI_z(ir, j, k)

    ir -> 0 for real, 1 for imaginary
    j -> jth x value
    k -> kth mu value

    ** defaults to ir = 0, j = 0, k = 'all' **

    
    steady.render()
    perturbed.render()

    Renders to screen steady/perturbed current values

    steady.getI()
    perturbed.getI()

    returns I / dI array

* OTHER USEFUL INFO

    - Rendering and plotting are independent of each
    other, meaning the program can be run, then plotted,
    then run again. Then rendering it will show an
    overplotted view of the two runs

    - Important functionality files are in ./src/
    
        src
        ├── App.py -> not implemented
        ├── Perturbed.py -> Perturbed class (run and plot implementation)
        ├── Steady.py -> Steady class (run and plot implementation)
        └── util
            ├── MatrixPerUtil.py -> Relax. method specific matrix creation
            ├── MatrixStUtil.py -> Relax. method specific matrix creation
            ├── PerUtil.py -> Perturbed specific source functions
            ├── StUtil.py -> Steady specific source functions
            ├── Util.py -> Shared source/other utility functions



********************************************** """

# EXAMPLE IMPLEMENTATION

steady = Steady.Steady()
steady.run(verb = 'v')
steady.plotI_z()
steady.render()

# perturbed = Perturbed.Perturbed()
# perturbed.run(verb = 'v')
# perturbed.plotI_z()
# perturbed.render()