
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

U = [8.0, 0.0, 8.0] 
An = [2, 8, 2]
Bn = [1.0, 0.0]


type = "Finite Well"


function test()
   #= Checks if the finite well gets at least 1 bound state =#

    energyEigenfunctions = getAllBoundStates(U, An, 1e-2)

    if length(energyEigenfunctions) > 0
        print("Finite well test passed with $(length(energyEigenfunctions)) bound states")
        grid, psi = psiSim(U, energyEigenfunctions[4], An, Bn, 0.01)
        plotWavefunction(grid, psi, energyEigenfunctions[4])
    else
        print("Finite well test failed")
    end
end

function tprpSim()
    energyEigenfunctions, tp, rp = getAllBoundStates(U, An, 1e-2)
    plotTprp(tp, rp)
end

function psiSim()
    energyEigenfunctions = getAllBoundStates(U, An, 1e-2)
    grid, psi = psiSim(U, energyEigenfunctions[4], An, Bn, 0.01)
    plotWavefunction(grid, psi, energyEigenfunctions[4])
end


function unboundStates()
end

test()