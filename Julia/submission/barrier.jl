

#= 
Bound States
Wavefunction for the bound states
Arbitary Potenital bound states
3 ev well of width 8 Angstrom
=#

include("constants.jl")
include("wavefunction.jl")
include("bound_states.jl")

energyEigenfunctions = []

U = [0, 2, 0] 
An = [2, 5, 2]
Bn = [1.0, 0.0]


type = "Barrier"


function test()
    #= Checks if the finite well gets at least 1 bound state =#
 
     energyEigenfunctions = getAllBoundStates(U, An, 1e-2)
 
     if length(energyEigenfunctions) > 0
        print("Barrier test passed with $(length(energyEigenfunctions)) bound states")
        grid, psi = psiSim(U, energyEigenfunctions[2], An, Bn, 0.01)
        plotWavefunction(grid, psi, energyEigenfunctions[2])
     else
        print("Barrier test failed")
     end
 end
 
function unboundStates()
end

test()