
function checkBoundStates(energyBoundStates, U, boundaries)
    """ Checks if boundstate array has pairs that reflect sign of t11 """
    
    i = 1
    ptr = true
    while i<length(energyBoundStates)
        
        lowerEnergy, higherEnergy = energyBoundStates[i], energyBoundStates[i+1]
        
        lowerEnergyt11, higherEnergyt11 = real(nReigonTm(lowerEnergy, U, boundaries)[1,1]), real(nReigonTm(higherEnergy, U, boundaries)[1,1])
        
        if lowerEnergyt11 <= 0.0 && higherEnergyt11 <= 0.0
            return false
        
        elseif lowerEnergyt11 >= 0.0 && higherEnergyt11 >= 0.0
            return false
        
        else
            i+=2
        end
    end

    return true
end

function boundStates(E, t11arr)
    #= Takes an input of the array of the t11 values and plots the t11 and returns the bound state energies =#

    # I think my midpoint function is fine, its just that t11 bound states does not 
    # correspond to a sign flip - need more redundency - definitely a problem here - write some tests for this
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
        
        # Getting midpoint from higher and lower energy
        midpoint = (lowerEnergy + higherEnergy)/2
        #println(lowerEnergy, "", higherEnergy, "", midpoint)
        
        t11 = nReigonTm(midpoint, U, boundaries)[1,1]
        println(t11)
        
        if abs(real(t11)) >= 1e-22
        
            lowerEnergyt11 = real(nReigonTm(lowerEnergy, U, boundaries)[1,1])
            higherEnergyt11 = real(nReigonTm(higherEnergy, U, boundaries)[1,1])
            
            # Finding which energy level has a larger t11 value and adjusting energy levels accordingly
            if real(t11) < 0  
                if lowerEnergyt11 > higherEnergyt11
                    higherEnergy = midpoint
                else
                    lowerEnergy = midpoint
                end
            else
                if lowerEnergyt11 > higherEnergyt11
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
    if checkBoundStates(energyBoundStates,U, boundaries)
    
        while n*2 <= length(energyBoundStates)
            lowerEnergy, higherEnergy = energyBoundStates[n], energyBoundStates[n+1]
            append!(energyEigenfunctions, getNthBoundStateEnergy(lowerEnergy, higherEnergy, U, boundaries))
            n = n + 1
        end
    return energyEigenfunctions

    else
        println("Bound states do not have same sign")
        return false
    end
end


