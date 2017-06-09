Relaxation method solver
========================
## Problem
Attempt to a solution for a perturbed linearized radiative transfer equation with coherent scattering terms.

## Method
Iterative solution that takes an initial guess for an intensity grid to solve for the difference between each point of the grid and the actual solution, and add this ∆I to the initial guess.

## Installation

Make sure a package manager such as `brew` or `apt-get` is available

If *Python 3* and/or *git* are **not** installed then install w/ dependencies using:
	
	brew install python
	brew install python3
	brew install git
	pip3 install numpy
	pip3 install matplotlib
	
> For *Linux* systems replace `brew install` with `sudo apt-get` 

- To download and install on *MacOS/Linux*:

		git clone https://github.com/pboserman/Relaxation.git
		cd Relaxation

- To allow updating:

	* Create a GitHub account.
	* On Repostory page click on **Fork** to fork to your profile
	* Run:

			git clone https://github.com/{username}/Relaxation.git 
			cd Relaxation
			git remote add upstream https://github.com/pboserman/Relaxation.git

	* Then to update local files:
	
			git pull upstream master
			
## Recommended installation

For interactive usage download *iPython*:

	brew install ipython

Then go to root directory and run:

	ipython
	run index.py
	
> *iPython* is an interactive *Python* shell that allows access to existing variables within the scope of the file ran.
> This is very useful for debugging and testing
	

## Usage
* To run the program, go to to root folder and run:
		
		python3 index.py

* Options can be set in the *opt.py* file in the root directory. 

* In order to make custom runs, look at the commented available methods within *index.html*.


## Other useful info

![Alt src] (images/src.jpg)	 *Source dir. information*