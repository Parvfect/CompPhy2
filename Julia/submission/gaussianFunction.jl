
include("potentials.jl")
include("bound_states.jl")

new_gaussian = createGaussianPotential(2.5, 0.2)
Bcn = [1.0, 0.0]
U, An = piecewisePotential(0, 5, 101, new_gaussian)


function test()
    #= Checks if the finite well gets at least 1 bound state =#
 
    energyEigenfunctions = getAllBoundStates(U, An, 1e-2)
    
    if length(energyEigenfunctions) > 0
        print("Gaussian test passed with $(length(energyEigenfunctions)) bound states")
        grid, psi = psiSim(U, energyEigenfunctions[2], An, Bn, 0.01)
        #plotWavefunction(grid, psi, energyEigenfunctions[2])
    else
        print("Gaussian test failed")
    end
 end

 test()