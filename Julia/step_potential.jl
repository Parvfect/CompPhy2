
include("wavefunction_1d.jl")

# Attempting to use my generic method to reproduce the results of the unit step funciton
e=1.6e-19
me=9.11e-31
# a=3e-10
hbar=1.05e-34

U = [0, 2] * e
E = [0.75, 1.5, 2.5] * e
Bc2 = [1.0, 0.0]
size_reigon = 3e-10 + 2e-9
An = [size_reigon, size_reigon]

grid, psi = nRegion(U, E[1], An, Bc2, 1e-11)

plotSimulation(grid, psi)