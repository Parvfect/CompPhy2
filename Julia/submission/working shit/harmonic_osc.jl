using Plots
using LinearAlgebra

e, me, hbar, A3 = 1.6e-19, 9.11e-31, 1.05e-34, 1.0
save_path = "C:/Users/Parv/Documents/compphy/Julia/Data/"
#U, reigon_lengths, boundaries = [3, 0, 3]*e, [2e-9, 8e-9, 2e-9], [2e-9, 10e-9]

function makePositive(x)
    if x < 0
        return -x
    else
        return x
    end
end

function getTransmissionProbability(t11, kn, k1)
    return makePositive(real((1/(t11))^2*(kn/k1)))
end

function getReflectionProbability(t21, t11)
    return makePositive(real((t21/t11)^2))
end

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

function nReigontprp(E, U, boundaries)
    TM, k = getTransferMatrix(getWaveVector(E, U[1]), getWaveVector(E, U[2]), boundaries[1]), getTotalWaveVector(E, U)
    return getTransmissionProbability(TM[1,1], k[length(k)], k[2]), getReflectionProbability(TM[2,1], TM[1,1])
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


function energyLoopt11(E, U, a, b, boundaries)
    return [nReigont11(i, U, boundaries) for i in E]
end

function energyLooptprp(E, U, a, b, boundaries)
    tp, rp = complex(zeros(length(E))), complex(zeros(length(E)))
    for i in 1:length(E)
        tp[i], rp[i] = nReigontprp(E[i], U, boundaries)
    end
    return tp, rp
end
function t11Sim()
    E = 0:1e-21:20e-19
    t11 = energyLoopt11(E, U, 1.0, 0, boundaries)
    E = 0:0.01:20
    t11 = [makePositive(real(i)) for i in t11]
    p1 = plot(E, real(t11), xlabel = "E", ylabel = "T11")
    display(p1)
end

function tprpSim()
    E = 0:1e-21:6e-19
    tp, rp = energyLooptprp(E, U, 1.0, 0, boundaries)
    E = 0:0.01:6
    p1 = plot(E, real(tp), xlabel = "E", ylabel = "Transmission Probability")
    p1 = plot!(E, real(rp), xlabel = "E", ylabel = "R11")
    display(p1)
end
    
function nReigonSim(E, U, boundaries, reigon_lengths)
    A, B = nReigon(E, U, 1.0, 0, boundaries)
    nReigonPlot(A, B, getTotalWaveVector(E, U), reigon_lengths, false)
end

function getBoundaries(reigon_lengths)
    ptr = 0
    boundaries = zeros(length(reigon_lengths) - 1)
    for i in 1:(length(reigon_lengths)-1)
        boundaries[i] = ptr + reigon_lengths[i]
        ptr = ptr + reigon_lengths[i]
    end
    return boundaries
end

function createGaussianPotential(V, a)
    #= V = 3 , a = 5 =#
    return x -> V * exp(- (x-a)^2/6.25)
end


function piecewisePotential(upperLimit, n_reigons, potentialFunction)
    lengthReigon = upperLimit / n_reigons
    U = zeros(n_reigons)
    ptr = lengthReigon/2

    for i in 1:n_reigons
        U[i] =  potentialFunction(ptr)
        ptr += lengthReigon
    end

    return U, [lengthReigon for i in 1:n_reigons]
end

p = createGaussianPotential(3,5)
U, an = piecewisePotential(10, 9, p)
U = U*e
boundaries = getBoundaries(an*1e-9)
t11Sim()
#tprpSim()
#nReigonSim(0.75*e, U, boundaries, an*1e-9)
# Bound States - Highest to Lowest
#E = 3.843359e-19
#E = 1.2678747928e-19
#E = 7.1410000239e-20
#E = 1.4122357142399997e-20
#E = 3.531310235e-21
#E = 0.8828720751e-21

