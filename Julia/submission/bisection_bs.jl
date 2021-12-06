# https://en.wikipedia.org/wiki/Bisection_method

#=
1. Calculate c midpoint of interval
2. Caclulate function value at midpoint
3. if Convergence is satisfactory, energyBoundStates
4. Examine sign of f(c) and replae either a or b with c 
=#

function boundStates(E, t11arr)
    #= Takes an input of the array of the t11 values and plots the t11 and returns the bound state energies =#

    #display(plot(E,t11arr))
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

    while true
        midpoint = (lowerEnergy + higherEnergy)/2
        t = getT11(midpoint, U, An)
        if abs(t) >= 0.001
            if t < 0 
                lowerEnergy = midpoint
            else
                higherEnergy = midpoint
            end
        else 
            println("Energy Eigenfunction found at $(E[i]) ev")
            return midpoint
        end
    end
end

function getAllBoundStates(U, E, An, acceptableError, step)

    energyEigenfunctions = []
    t11arr, k_arr, A_arr, B_arr = energyLoop(E, U, An, step)
    n = 1
    println(A_arr)
    energyBoundStates = boundStates(E, t11arr)
    while n*2 <= length(energyBoundStates)
        append!(energyEigenfunctions, getNthBoundStateEnergy(n, E, t11arr, U, An, energyBoundStates,acceptableError, step))
        n = n + 1
    end
    return energyEigenfunctions
end


