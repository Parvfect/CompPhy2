
include("potentials.jl")
include("bound_states.jl")
include("useful_methods.jl")

new_gaussian = createGaussianPotential(2, 5)
Bcn = [1.0, 0.0]
U, An = piecewisePotential(10, 101, new_gaussian)
step = 1e-4
E = createRange(0.01, 2.99, 1e-2)

function test()
    #= Checks if the finite well gets at least 1 bound state =#

    energyEigenfunctions = getAllBoundStates(U, E, An, 1e-2, 1e-2)

    if length(energyEigenfunctions) > 0
        print("Gaussian test passed with $(length(energyEigenfunctions)) bound states")
        A, B, K,
        t = solveTMM(U, energyEigenfunctions[2], An, Bn, 1e-2)
        psi = totalWavefunction(An, A, B, K)
        #grid = 0:0.01:12.02
        #plotWavefunction(grid, psi, 1.13)
    else
    print("Gaussian test failed")
    end
end
 
test()