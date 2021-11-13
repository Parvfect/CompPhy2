
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

parameters = [
    U = [3, 0, 3] 
    An = [4, 8, 4]
    Bn = [1.0, 0.0]
]

type = "Finite Well"


function test()
   #= Checks if the finite well gets at least 1 bound state =#

    energyEigenfuncitons = getAllBoundStates(parameters)

    if length(energyEigenfuncitons) > 0
        print("Finite well test passed with $(length(energyEigenfunctions)) bound states")
        plotWavefunctions(energyEigenfunctions, type)
    else
        print("Finite well test failed")
    end
end

function unboundStates()
end