
include("wavefunction_1d.jl")

# Attempting to use my generic method to reproduce the results of the unit step funciton
e=1.6e-19
me=9.11e-31
# a=3e-10
hbar=1.05e-34

U = [0, 2] * e
E = [0.75, 1.5, 2.5] * e
Bc2 = [1.0, 0.0]
An = [23e-11, 3 + 20]

grid, psi = nRegion(U, E[1], An, Bc2, -20, 20, 0.01)
#plotSimulation(grid, psi)
grid = reverse(grid)
psi = reverse(psi)
plotSimulation(psi, grid)