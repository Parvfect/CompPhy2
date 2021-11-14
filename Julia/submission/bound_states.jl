
include("useful_methods.jl")
include("wavefunction.jl")
include("constants.jl")
using Base

function boundStates(E, t11arr)
    #= Takes an input of the array of the t11 values and plots the t11 and returns the bound state energies =#

    # Plot the t11 values over E
    display(plot!(E, t11arr, xlabel = "Energies(eV) ", ylabel = "t11", title = "Variation of the t11 over Energy Range"))
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

function getNthBoundStateEnergy(n, E, t11arr, U, An, acceptableError)
    #= Gets the nth bound state Sign crossover for t11. Repeats simulation between the two energies till energy eigenfunction is found =#
        
    energyBoundStates = boundStates(E, t11arr)
    if n*2 > length(energyBoundStates)
        println("Error: n is too large")
        return -1
    end

    # Get the energy of the nth bound state
    lowerEnergy = energyBoundStates[n*2-1]
    higherEnergy = energyBoundStates[n*2]

    if abs(t11arr[lowerEnergy]) <= acceptableError
        println("Energy Eigenfunction found at $lowerEnergy ev")
        return lowerEnergy

    elseif abs(t11arr[higherEnergy]) <= acceptableError
        println("Energy Eigenfunction found at $higherEnergy ev")
        return higherEnergy
    end
        

    iterations = 0

    # Find Energy Eigenfunction between the two energies
    while iterations < 10
    
        E = createRange(lowerEnergy, higherEnergy, step/(10^iterations))

        t11 = energyLoop(E, U, An)

        for i in 1:length(t11)
            if abs(t11[i]) <= accetableError 
                println("Energy Eigenfunction found at $(E[i]) ev")
                return E[i]
            end
        end
    end
    
    println("Error: Could not find the energy of the nth bound state")
    return -1
end

function energyLoop(E, U, An)
    
    t11arr = zeros(0)

    @time for i in E
        t11arr = append!(t11arr, real(t11Sim(U, i, An, [1.0, 0.0], 1e-11)))
    end

    return t11arr
end


function getAllBoundStates(U, An, acceptableError)

    energyEigenfunctions = []

    E = createRange(0.0, maximum(U), 1e-11)

    t11arr = energyLoop(E, U, An)

    n = 1

    while n*2 <= length(boundStates(E, t11arr))
        append!(energyEigenfunctions, getNthBoundStateEnergy(n, E, t11arr, acceptableError, U, An))
        n = n + 1
    end

    for i in 1:length(energyEigenfunctions)
        if energyEigenfunctions[i] == -1
            println("Error: Could not find the energy of the $i th bound state")
            return -1
        end
    end

    return energyEigenfunctions
end

function boundStatesCheck()
    #= Checks the bound states of the t11 values =#

    t11 = [2,1,0,-1,-2,-1,0,1,2,1,0,-1,-2.3e-5,0]
    E = [0,1,2,3,4,5,6,7,8,9,10,11,12,13]

    en = boundStates(E, t11)

    if en == [2, 3, 5, 6, 10, 11,12,13]
        println("bound States finder test passed")
    else
        println("bound States finder test failed")
    end
end

function getNthBoundStateEnergyCheck()
    #= Checks the getNthBoundStateEnergy function =#

    t11 = [2,1,0,-1,-2,-1,0,1,2,1,0,-1,-2.3e-5,0]
    E = [0,1,2,3,4,5,6,7,8,9,10,11,12,13]

    en = boundStates(E, t11)

    if getNthBoundStateEnergy(1, E, t11, 1e-3, U, An) == 3
        println("getNthBoundStateEnergy test passed")
    else
        println("getNthBoundStateEnergy test failed")
    end
end
