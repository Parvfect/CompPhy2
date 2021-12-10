using Plots
using LinearAlgebra

e, me, hbar, A3 = 1.6e-19, 9.11e-31, 1.05e-34, 1.0
save_path = "C:/Users/Parv/Documents/compphy/Julia/Data/"
U, reigon_lengths, boundaries = [3, 0, 3]*e, [2e-9, 8e-9, 2e-9], [2e-9, 10e-9]


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

function nReigont11(E, U, boundaries)
    return getTransferMatrix(getWaveVector(E, U[1]), getWaveVector(E, U[2]), boundaries[1])[1,1]       
end

function nReigon(E, U, a, b, boundaries)

    A, B, k = complex(zeros(length(U))), complex(zeros(length(U))), getTotalWaveVector(E, U)
    A[length(U)], B[length(U)] = a, b

    for i in length(A):-1:2
        TM = getTransferMatrix(k[i-1], k[i], boundaries[i-1])
        AB = TM*[A[i];B[i]]
        A[i-1], B[i-1] = AB[1], AB[2]
    end
    return A, B
end

function nReigonPlot(A, B, k, reigon_lengths, saveFig)
    
    ptr = 0
    p1 = plot(xlabel = "x", ylabel = "psi(x)", xlim=(0,1.2e-8), ylim=(-Inf,+Inf), grid = true)
    for i in 1:length(A)
        x = ptr:1e-11:ptr + reigon_lengths[i]
        psi = [getWavefunction(A[i], B[i], k[i], j) for j in x]
        p1 = plot!(x, real(psi))
        ptr = ptr + reigon_lengths[i]
    end
    display(p1)

    if saveFig
        s = readline("Enter Figure name")
        savefig(p1, save_path*"$s")
    end
end


function energyLoop(E, U, a, b, boundaries)
    return [nReigont11(i, U, boundaries) for i in E]
end

function t11Sim()
    E = 0:1e-21:3*e
    t11 = energyLoop(E, U, 1.0, 0, boundaries)
    p1 = plot(real(t11), xlabel = "E", ylabel = "T11")
    display(p1)
end

function nReigonSim(E)
    A, B = nReigon(E, U, 1.0, 0, boundaries)
    nReigonPlot(A, B, getTotalWaveVector(E, U), reigon_lengths, false)
end


# Bound States - Highest to Lowest
#E = 3.843359e-19
#E = 1.2678747928e-19
#E = 7.1410000239e-20
E = 1.4122357142399997e-20
#E = 3.531310235e-21
#E = 0.8828720751e-21

t11Sim()