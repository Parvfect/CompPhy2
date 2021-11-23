
#= 
Bound States
Wavefunction for the bound states
Arbitary Potenital bound states
3 ev well of width 8 Angstrom
=#

include("constants.jl")
include("wavefunction.jl")
include("bound_states.jl")
include("useful_methods.jl")

energyEigenfunctions = []

U = [3.0, 0.0, 3.0] 
An = [2, 8, 2]
Bn = [1.0, 0.0]
E = createRange(0.01, 2.99, 1e-2)

type = "Finite Well"

function test()
   #= Checks if the finite well gets at least 1 bound state =#

    energyEigenfunctions = getAllBoundStates(U, E, An, 1e-2, 1e-2)

    if length(energyEigenfunctions) > 0
        print("Finite well test passed with $(length(energyEigenfunctions)) bound states")
            A, B, K, t = solveTMM(U, energyEigenfunctions[2], An, Bn, 1e-2)
            psi = totalWavefunction(An, A, B, K)
            #grid = 0:0.01:12.02
            #plotWavefunction(grid, psi, 1.13)
    else
        print("Finite well test failed")
    end
end

function tprpSim()
    energyEigenfunctions = getAllBoundStates(U, An, 1e-2)
end

function psiSim()
    grid, psi = psiSim(U, 1.08999, An, Bn, 0.01)
    plotWavefunction(grid, psi, 1.35)
end


function unboundStates()
end

test()