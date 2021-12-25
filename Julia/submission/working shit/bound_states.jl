#=  .../bound_states.jl 
    Finds the energy bound States of any system given 
    a function that returns the t11 for a certain energy
    and a t11 that intersects zero
=#


function checkFlip(lowerEnergy, higherEnergy, U, boundaries)
    """ Checks if boundstate array has pairs that reflect sign of t11 """
    
    lowerEnergyt11, higherEnergyt11 = real(nReigonTm(lowerEnergy, U, boundaries)[1,1]), real(nReigonTm(higherEnergy, U, boundaries)[1,1])
    
    if lowerEnergyt11 <= 0.0 && higherEnergyt11 <= 0.0 || lowerEnergyt11 >= 0.0 && higherEnergyt11 >= 0.0
        return false
    else
        return true
    end
end

function boundStates(E, t11arr)
    #= Takes an input of the array of the t11 values and plots the t11 and returns the bound state energies =#

    oldValue = t11arr[1]
    energyBoundStates = zeros(0)    
    for i in 1:length(t11arr)
        if t11arr[i] >= 0 && oldValue < 0 || oldValue >= 0 && t11arr[i] < 0
            append!(energyBoundStates, E[i], E[i-1])
        end
        oldValue = t11arr[i]
    end
    return energyBoundStates
end

function getNthBoundStateEnergy(lowerEnergy, higherEnergy,  U, boundaries)
    #= Gets the nth bound state Sign crossover for t11. Repeats simulation between the two energies till energy eigenfunction is found =#
        
    while true
        
        if !checkFlip(lowerEnergy, higherEnergy, U, boundaries)
            throw(DomainError(ch, "Sign Flip lost - suggests that the function manipulated too much or there are no bound states"))
        end
        # Getting midpoint from higher and lower energy
        midpoint = (lowerEnergy + higherEnergy)/2
        #println(lowerEnergy, "", higherEnergy, "", midpoint)
        
        t11 = real(nReigonTm(midpoint, U, boundaries)[1,1])

        # Minimum error seems to be 1e-13 and then the thing starts repeating
        if abs(real(t11)) >= 1e-13
        
            lowerEnergyt11 = real(nReigonTm(lowerEnergy, U, boundaries)[1,1])
            higherEnergyt11 = real(nReigonTm(higherEnergy, U, boundaries)[1,1])
            
            # Finding which energy level has a larger t11 value and adjusting energy levels accordingly
            if higherEnergyt11 > lowerEnergyt11   
                # Since we know that there is a 0 between he and le we can use this feature to find the energy level
                if t11 > 0
                    higherEnergy = midpoint
                else
                    lowerEnergy = midpoint
                end
            else
                if t11 > 0
                    lowerEnergy = midpoint
                else
                    higherEnergy = midpoint
                end
            end

        # If t11 is less than absolute error, return the bound state   
        else 
            println("Energy Eigenfunction found at $midpoint ev")
            return midpoint
        end
    end
end

function getAllBoundStates(E, U, boundaries)

    energyEigenfunctions = []
    t11arr = real(energyLoop(E, U, boundaries))
    n = 1
    energyBoundStates = boundStates(E, t11arr)
    while n*2 <= length(energyBoundStates)
        lowerEnergy, higherEnergy = energyBoundStates[n], energyBoundStates[n+1]
        append!(energyEigenfunctions, getNthBoundStateEnergy(lowerEnergy, higherEnergy, U, boundaries))
        n = n + 1
    end
    return energyEigenfunctions
end


