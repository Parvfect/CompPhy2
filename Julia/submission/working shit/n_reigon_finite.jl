
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

function nReigon(E, U, a, b, boundaries)

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
    return A, B
end

function nReigonPlot(A, B, k, reigon_lengths)
    
    ptr = 0
    p1 = plot()
    for i in 1:length(A)
        x = ptr:1e-11:ptr + reigon_lengths[i]
        psi = [getWavefunction(A[i], B[i], k[i], j) for j in x]
        p1 = plot!(x, real(psi))
        ptr = ptr + reigon_lengths[i]
    end
    display(p1)
end


function t11plot(E)
    t11arr = getT11(E)
    display(plot(E,t11arr))
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


function t11plot(E)
    t11arr = getT11(E)
    display(plot(E,t11arr))
end


e=1.6e-19

E = 0:1e-21:2.99*e
U = [2, 0, 2]*e
reigon_lengths = [2e-9, 5e-9, 2e-9]
boundaries = [2e-9, 7e-9]
me=9.11e-31 
hbar=1.05e-34
A3=1.0
E = 1.318366*e
#A, B = nReigon(E, U, 1.0, 0, boundaries)

#println(solve(E))
A, B = nReigon(E, U, 1.0, 0, boundaries)
nReigonPlot(A, B, getTotalWaveVector(E, U), reigon_lengths)
# -559.4067088205346 - 4133.246567573868im - t11 for that bound state