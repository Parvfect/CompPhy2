

#= 
Let's just think here right
Finite well has three reigons 
We just need bound states for the finite well
Let's fuck the nreigon for now let's just do two reigons
manually and couple them together, and see if our bound states function works 
and gives us eigenfunctions =#

#= mark's solution for potential step
 
using Plots
using LinearAlgebra

e=1.6e-19

U1=0.0*e
U2=2*e

E=0.75*e

me=9.11e-31
a=3e-10
hbar=1.05e-34

k1=sqrt(2*me*(Complex(E-U1)))/hbar
k2=sqrt(2*me*(Complex(E-U2)))/hbar

A2=1.0

TM=(1/(2*k1))*[(k1+k2)*exp(-1im*a*(k1-k2)) (k1-k2)*exp(-1im*a*(k1+k2)); (k1-k2)*exp(1im*a*(k1+k2))  (k1+k2)*exp(1im*a*(k1-k2))]*[A2;0]
A1=TM[1]
B1=TM[2]

x=-2e-9:1e-11:a

psi1=A1*exp.(1im.*x.*k1).+B1*exp.(-1im.*x.*k1)

p1=plot(x,real(psi1))
p1=plot!(x,imag(psi1))

x=a:1e-11:2e-9

psi2=A2*exp.(1im.*x.*k2)

p1=plot!(x,real(psi2))
p1=plot!(x,imag(psi2))

=# 

#= Generalising Mark's solution for a x>0 solution 
using Plots
using LinearAlgebra

e=1.6e-19

U1=0.0*e
U2=2*e

E=0.75*e

me=9.11e-31
a=5e-9
hbar=1.05e-34

k1=sqrt(2*me*(Complex(E-U1)))/hbar
k2=sqrt(2*me*(Complex(E-U2)))/hbar

A2=1.0

TM=(1/(2*k1))*[(k1+k2)*exp(-1im*a*(k1-k2)) (k1-k2)*exp(-1im*a*(k1+k2)); (k1-k2)*exp(1im*a*(k1+k2))  (k1+k2)*exp(1im*a*(k1-k2))]*[A2;0]
A1=TM[1]
B1=TM[2]

x=0:1e-11:a

psi1=A1*exp.(1im.*x.*k1).+B1*exp.(-1im.*x.*k1)

p1=plot(x,real(psi1))
p1=plot!(x,imag(psi1))

x=a:1e-11:10e-9

psi2=A2*exp.(1im.*x.*k2)

p1=plot!(x,real(psi2))
p1=plot!(x,imag(psi2))

=#

#= Making Potential Step work with functions 
function getWaveVector(E, U)
    return sqrt(2*me*(Complex(E-U)))/hbar
end

function getTransferMatrix(k1, k2, a)
    return (1/(2*k1))*[(k1+k2)*exp(-1im*a*(k1-k2)) (k1-k2)*exp(-1im*a*(k1+k2)); (k1-k2)*exp(1im*a*(k1+k2))  (k1+k2)*exp(1im*a*(k1-k2))]
end

function getWavefunction(A, B, k, x)
    return A*exp.(1im.*x.*k).+B*exp.(-1im.*x.*k)
end

using Plots
using LinearAlgebra

e=1.6e-19

U1=0.0*e
U2=2*e

E=0.75*e

me=9.11e-31
a1=2e-9
a2=10e-9
hbar=1.05e-34

k1=getWaveVector(E, U1)
k2=getWaveVector(E, U2)

A2=1.0

TM=getTransferMatrix(k1, k2, a)
AB = TM*[A2;0]
A1=TM[1]
B1=TM[2]

x=0:1e-11:a

psi1=getWavefunction(A1, B1, k1, x)

p1=plot(x,real(psi1))
p1=plot!(x,imag(psi1))

x=a:1e-11:10e-9

psi2=getWavefunction(A2, 0, k2, x)

p1=plot!(x,real(psi2))
p1=plot!(x,imag(psi2))

=#
#= Making functions out of Mark's solution 

using Plots
using LinearAlgebra

e=1.6e-19

U1=2*e
U2=0*e
U3 = 2*e

me=9.11e-31
a1=2e-9
a2 = 10e-9 
hbar=1.05e-34
A3 = 1.0

function getWaveVector(E, U)
    return sqrt(2*me*(Complex(E-U)))/hbar
end

function getTransferMatrix(k1, k2, a)
    return (1/(2*k1))*[(k1+k2)*exp(-1im*a*(k1-k2)) (k1-k2)*exp(-1im*a*(k1+k2)); (k1-k2)*exp(1im*a*(k1+k2))  (k1+k2)*exp(1im*a*(k1-k2))]
end

function getLeftBoundary(k1, k2, A, B, a)
    return (1/(2*k1))*[(k1+k2)*exp(-1im*a*(k1-k2)) (k1-k2)*exp(-1im*a*(k1+k2)); (k1-k2)*exp(1im*a*(k1+k2))  (k1+k2)*exp(1im*a*(k1-k2))]*[A;B]
end

function getWavefunction(A, B, k, x)
    return A*exp.(1im.*x.*k).+B*exp.(-1im.*x.*k)
end

function getT11(E)
    k1 = getWaveVector(E, U1)
    k2= getWaveVector(E, U2)
    k3 = getWaveVector(E, U3)
    TM=getLeftBoundary(k2, k3, A3, 0, a2)
    #A2=TM[1]
    #B2=TM[2]
    #TM=getTransferMatrix(k1, k2, a1)
    #t11 = TM[1,1]
    return TM
end

function solve(E)

    E = E*e
    k1 = getWaveVector(E, U1)
    k2= getWaveVector(E, U2)
    k3 = getWaveVector(E, U3)
    
    TM=getLeftBoundary(k2, k3, A3, 0, a2)
    A2=TM[1]
    B2=TM[2]
    
    
    TM=getLeftBoundary(k1, k2, A2, B2, a1)
    A1=TM[1]
    B1=TM[2]
    t11 = TM[1,1]
    
    
    x1=0:1e-11:a1
    
    psi1=getWavefunction(A1, B1, k1, x1)
    
    p1=plot(x1,real(psi1))
    #p1=plot!(x,imag(psi1))
    
    x2=a1:1e-11:a2
    
    psi2=getWavefunction(A2, B2, k2, x2)
    
    p1=plot!(x2,real(psi2))
    #p1=plot!(x,imag(psi2))
     
    x3=a2:1e-11:12e-9
    
    psi3=getWavefunction(A3, 0, k3, x3)
    p1=plot!(x3,real(psi3))
    #p1=plot!(x,imag(psi3))

    display(p1)
end

solve(0.75)
=#
#=
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
    t11 = [getT11[i] for i in E]
    n = 1
    println(A_arr)
    energyBoundStates = boundStates(E, t11arr)
    while n*2 <= length(energyBoundStates)
        append!(energyEigenfunctions, getNthBoundStateEnergy(n, E, t11arr, U, An, energyBoundStates,acceptableError, step))
        n = n + 1
    end
    return energyEigenfunctions
end

=#

#= Finite well through new solution method (without generalising) 
function getWaveVector(E, U)
    return sqrt(2*me*(Complex(E-U)))/hbar
end

function getTransferMatrix(k1, k2, a)
    return (1/(2*k1))*[(k1+k2)*exp(-1im*a*(k1-k2)) (k1-k2)*exp(-1im*a*(k1+k2)); (k1-k2)*exp(1im*a*(k1+k2))  (k1+k2)*exp(1im*a*(k1-k2))]
end

function getWavefunction(A, B, k, x)
    return A*exp.(1im.*x.*k).+B*exp.(-1im.*x.*k)
end

using Plots
using LinearAlgebra

e=1.6e-19

U1=3.0*e
U2=0*e
U3=3*e
E=1.025*e

me=9.11e-31
a1=2e-9
a2=10e-9
hbar=1.05e-34

k1=getWaveVector(E, U1)
k2=getWaveVector(E, U2)
k3=getWaveVector(E, U3)

A3=1.0

TM=getTransferMatrix(k2, k3, a2)
AB = TM*[A2;0]
A2=TM[1]
B2=TM[2]

TM=getTransferMatrix(k1, k2, a1)
print(TM)
AB = TM*[A2;0]
A1=TM[1]
B1=TM[2]


x=0:1e-11:a1

psi1=getWavefunction(A1, B1, k1, x)

p1=plot(x,real(psi1))
#p1=plot!(x,imag(psi1))

x=a1:1e-11:a2

psi2=getWavefunction(A2, B2, k2, x)

p1=plot!(x,real(psi2))
#p1=plot!(x,imag(psi2))

x=a2:1e-11:12e-9
psi3 = getWavefunction(A3, 0, k3, x)
p1=plot!(x,real(psi3))
#p1=plot!(x,imag(psi3))
=#

#= Generalising Finite well

using Plots
using LinearAlgebra


e=1.6e-19

U = [3.0, 0, 3.0]*e
E=1.025*e
reigon_lengths = [2e-9, 8e-9, 2e-9]
boundaries = [2e-9, 10e-9]
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
    A2=TM[1]
    B2=TM[2]

    TM=getTransferMatrix(k1, k2, boundaries[1])
    AB = TM*[A2;B2]
    A1=TM[1]
    B1=TM[2]

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
end

E = 0:1e-21:2*e
t11 = getT11(E)
display(plot(t11))
=#

#= Bound States of Finite Well  - yum


using Plots
using LinearAlgebra


e=1.6e-19

U = [0, 2, 0]*e
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
solve(1.1*e)

#Looks like we are getting the right solution for 
#the potential barrier - who knows what the fuck 


=#

# Generalising through n reigon simulation (let's get this bread) 

using Plots
using LinearAlgebra

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

function getTotalWaveVector(E, U)
    return [getWaveVector(E, i) for i in U]
end

function energyLoop(E, U, a, b, boundaries)
    return [nReigont11(i, U, a, b, boundaries) for i in E]
end

function nReigont11(E, U, a, b, boundaries)
    k = getTotalWaveVector(E, U)
    A = complex(zeros(length(k)))
    B = complex(zeros(length(k)))
    A[length(k)] = a
    B[length(k)] = b
    t11 = 0 

    for i in length(A):-1:2
    
        TM = getTransferMatrix(k[i-1], k[i], boundaries[i-1])
        if i == 2
            t11 = TM[1,1]
        end
        AB = TM*[A[i];B[i]]
        A[i-1] = AB[1] 
        B[i-1] = AB[2]
    end
    return t11
end

function nReigon(E, U, a, b, boundaries)

    k = getTotalWaveVector(E, U)
    A = complex(zeros(length(k)))
    B = complex(zeros(length(k)))
    A[length(k)] = a
    B[length(k)] = b

    for i in length(A):-1:2
    
        TM = getTransferMatrix(k[i-1], k[i], boundaries[i-1])
        AB = TM*[A[i];B[i]]
        A[i-1] = AB[1] 
        B[i-1] = AB[2]
    end
    return A, B
end

function nReigonPlot(A, B, k, reigon_lengths)
    
    ptr = 0
    p1 = plot(xlabel = "x", ylabel = "psi(x)", xlims = 0)
    for i in 1:length(A)
        x = ptr:1e-11:ptr + reigon_lengths[i]
        psi = [getWavefunction(A[i], B[i], k[i], j) for j in x]
        println(x)
        p1 = plot!(x, real(psi))
        
        ptr = ptr + reigon_lengths[i]
    end
    display(p1)
    
    #savefig(p1,"C:/Users/Parv/Documents/compphy/Julia/Data/BarrierWavefunction")
end

function t11plot(E, U, boundaries)
    t11 = energyLoop(E, U, 1.0, 0, boundaries)
    p1 = plot(xlabel = "E", ylabel = "T11")    
    p1 = plot!(real(t11))
    display(p1)
    savefig(p1,"C:/Users/Parv/Documents/compphy/Julia/Data/t11arrFW")
end

function boundStates(E, t11arr)
    #= Takes an input of the array of the t11 values and plots the t11 and returns the bound state energies =#

    #display(plot(E,t11arr))
    println(E)
    println(t11arr)
    t11arr = real(t11arr)
    oldValue = t11arr[1]
    energyBoundStates = zeros(0)    
    for i in 1:length(t11arr)
        if t11arr[i] >= 0 && oldValue < 0 || oldValue >= 0 && t11arr[i] < 0
            append!(energyBoundStates, E[i-2], E[i-1])
            #println(oldValue, t11arr[i])
            #println(E[i-1], E[i])
        end
        oldValue = t11arr[i]
    end
    return energyBoundStates
end

function getNthBoundStateEnergy(n, E, t11arr, U, A, B, energyBoundStates, boundaries)
    #= Gets the nth bound state Sign crossover for t11. Repeats simulation between the two energies till energy eigenfunction is found =#
        
    lowerEnergy = energyBoundStates[n*2-1]
    higherEnergy = energyBoundStates[n*2]    
    if abs(real(nReigont11(lowerEnergy, U, A, B, boundaries))) < 0.01
        return lowerEnergy
    end
    if abs(real(nReigont11(higherEnergy, U, A, B, boundaries))) < 0.01
        return higherEnergy
    end
    while true
        midpoint = (lowerEnergy + higherEnergy)/2 
        t = real(nReigont11(midpoint, U, A, B, boundaries))
        if abs(t) > 0.01
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

function getAllBoundStates(U, E, A, B, boundaries)

    energyEigenfunctions = []
    t11arr = real(energyLoop(E, U, A, B, boundaries))
    display(plot(E, real(t11arr)))
    n = 1
    energyBoundStates = boundStates(E, t11arr)
    println(energyBoundStates)
    while n*2 <= length(energyBoundStates)
        append!(energyEigenfunctions, getNthBoundStateEnergy(n, E, t11arr, U, A, B, energyBoundStates, boundaries))
        n = n + 1
    end
    return energyEigenfunctions
end


e=1.6e-19
save_path = "C:/Users/Parv/Documents/compphy/Julia/Data/"
U = [3, 0, 3]*e
reigon_lengths = [2e-9, 8e-9, 2e-9]
boundaries = [2e-9, 10e-9]
me=9.11e-31 
hbar=1.05e-34
A3=1.0
#E = 2.4021*e
#A, B = nReigon(E, U, 1.0, 0, boundaries)
E = 0:1e-21:2.9e-19
println(length(E))
t11 = energyLoop(E, U, 1.0, 0, boundaries)
#for i in 1:(length(E))
 #   println(E[i], t11[i])
#end
#E = 0:0.01:2.90
#println(length(E))
#display(plot(E, real(t11)))
#println(solve(E))
#A, B = nReigon(E, U, 1.0, 0, boundaries)

#new_gaussian = createGaussianPotential(2, 5)
#U, reigon_lengths = piecewisePotential(10, 101, new_gaussian)
#boundaries = []
#ptr = 0

#=
for i in 1:(length(reigon_lengths)-1)
    boundaries = append!(boundaries, ptr + reigon_lengths[i])
    ptr += reigon_lengths[i]
end
print(boundaries)
t11plot(E, U, boundaries)
=#
#nReigonPlot(A, B, getTotalWaveVector(E, U), reigon_lengths)
# -559.4067088205346 - 4133.246567573868im - t11 for that bound state

print(getAllBoundStates(U, E, 1.0, 0, boundaries))