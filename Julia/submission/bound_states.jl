
include("useful_methods.jl")
include("wavefunction.jl")
include("constants.jl")
using Base

function boundStates(E, t11arr)
    #= Takes an input of the array of the t11 values and plots the t11 and returns the bound state energies =#

    oldValue = t11arr[1]
    energyBoundStates = zeros(0)    
    for i in 1:length(t11arr)
        if t11arr[i] >= 0 && oldValue < 0 || oldValue >= 0 && t11arr[i] < 0
            append!(energyBoundStates, E[i-1], E[i])
        end
        oldValue = t11arr[i]
    end
    return energyBoundStates
end

function getNthBoundStateEnergy(n, E, t11arr, U, An, energyBoundStates, acceptableError, step)
    #= Gets the nth bound state Sign crossover for t11. Repeats simulation between the two energies till energy eigenfunction is found =#
        
    lowerEnergy = energyBoundStates[n*2-1]
    higherEnergy = energyBoundStates[n*2]    

    iterations = 0
    while iterations < 10
        E = createRange(lowerEnergy, higherEnergy, 0.01/(10^iterations))
        t11, k_arr, A_arr, B_arr = energyLoop(E, U, An, step)
        for i in 1:length(t11)
            if abs(t11[i]) <= acceptableError 
                println("Energy Eigenfunction found at $(E[i]) ev")
                return E[i]
            end
        end
        iterations += 1
    end
    println("Error: Could not find the energy of the nth bound state")
    return -1
end

function getAllBoundStates(U, E, An, acceptableError, step)

    energyEigenfunctions = []
    t11arr, k_arr, A_arr, B_arr = energyLoop(E, U, An, step)
    n = 1
    energyBoundStates = boundStates(E, t11arr)
    while n*2 <= length(energyBoundStates)
        append!(energyEigenfunctions, getNthBoundStateEnergy(n, E, t11arr, U, An, energyBoundStates,acceptableError, step))
        n = n + 1
    end
    return energyEigenfunctions
end

