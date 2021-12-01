
using Plots
using LinearAlgebra


e=1.6e-19

U = [2, 0, 2]*e
E=1.025*e
reigon_lengths = [2e-9, 5e-9, 2e-9]
boundaries = [2e-9, 7e-9]
me=9.11e-31 
hbar=1.05e-34
A3=1.0

function getWaveVector(E, U)
    return sqrt(2*me*(Complex(E-U)))/hbar
end

function getTransferMatrix(k1, k2, a)
    return (1/(2*k1))*[(k1+k2)*exp(-1im*a*(k1-k2)) (k1-k2)*exp(-1im*a*(k1+k2)); (k1-k2)*exp(1im*a*(k1+k2))  (k1+k2)*exp(1im*a*(k1-k2))]
end

function getWavefunction(A, B, k, x)
    return A*exp.(1im.*x.*k).+B*exp.(-1im.*x.*k)
end

function getT11(E)
    return [real(solve(i)) for i in E]
end

function solve(E)
    k1=getWaveVector(E, U[1])
    k2=getWaveVector(E, U[2])
    k3=getWaveVector(E, U[3])


    TM=getTransferMatrix(k2, k3, boundaries[2])
    AB = TM*[A3;0]
    A2=AB[1]
    B2=AB[2]

    TM=getTransferMatrix(k1, k2, boundaries[1])
    AB = TM*[A2;B2]
    A1=AB[1]
    B1=AB[2]
    plot_simulation(A1, B1, A2, B2, A3, 0, k1, k2, k3)

    return TM[1,1]
end

function plot_simulation(A1, B1, A2, B2, A3, B3, k1, k2, k3)

    x=0:1e-11:boundaries[1]
    psi1=getWavefunction(A1, B1, k1, x)

    p1=plot(x,real(psi1))
    #p1=plot!(x,imag(psi1))

    x=boundaries[1]:1e-11:boundaries[2]

    psi2=getWavefunction(A2, B2, k2, x)

    p1=plot!(x,real(psi2))
    #p1=plot!(x,imag(psi2))

    x=boundaries[2]:1e-11:boundaries[2]+reigon_lengths[length(reigon_lengths)]
    psi3 = getWavefunction(A3, 0, k3, x)
    p1=plot!(x,real(psi3))
    #p1=plot!(x,imag(psi3))
    display(p1)
end

function createEnergyRange(lowerBound, stepDenom, upperBound)
    return lowerBound:(1e-21/10^stepDenom):upperBound
end

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

function getNthBoundStateEnergy(n, energyBoundStates, acceptableError)
    #= Gets the nth bound state Sign crossover for t11. Repeats simulation between the two energies till energy eigenfunction is found =#
        
    lowerEnergy = energyBoundStates[n*2-1]
    higherEnergy = energyBoundStates[n*2]    
    iterations = 0
    while iterations < 10
        E = createEnergyRange(lowerEnergy, iterations, higherEnergy)
        t11 = getT11(E)
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

function getAllBoundStates(E, acceptableError)

    energyEigenfunctions = []
    t11arr = getT11(E)
    n = 1
    
    energyBoundStates = boundStates(E, t11arr)
    while n*2 <= length(energyBoundStates)
        append!(energyEigenfunctions, getNthBoundStateEnergy(n, energyBoundStates, acceptableError))
        n = n + 1
    end
    return energyEigenfunctions
end


function t11plot(E)
    t11arr = getT11(E)
    display(plot(E,t11arr))
end

E = 0:1e-21:2.99*e
#t11 = [abs(i) for i in (getT11(E))]
#display(plot(t11)) 

#print(getAllBoundStates(E, 1e-3))
solve(1.318366*e)

#Looks like we are getting the right solution for 
#the potential barrier - who knows what the fuck 

