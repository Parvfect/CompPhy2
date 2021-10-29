
include("wavefunction_1d.jl")

# Attempting to use my generic method to reproduce the results of the unit step funciton

U = [0, 2] * e
E = [0.75, 1.5, 2.5] * e
Bc2 = [1.0, 0.0]
An = [3 + 20, 3 + 20]

grid, psi = nRegion(U, E[1], An, Bc2, -20, 20, 0.01)
#plotSimulation(grid, psi)
plotSimulation(psi[1:2000], grid[1:2000])