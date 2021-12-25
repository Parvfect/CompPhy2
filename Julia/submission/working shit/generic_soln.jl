using Plots
using LinearAlgebra

save_path = "C:/Users/Parv/Documents/compphy/Julia/Data/"

""" Helper Methods """

function getWaveVector(E, U)
    return sqrt(2*me*(Complex(E-U)))/hbar
end

function modulus(z)
    return sqrt(real(z) * real(z) + imag(z) * imag(z));
end

function getTransferMatrix(k1, k2, a)
    if k1 == 0
        print("k1 is zero")
        return [(1 - k2*a)*exp(1im *k2 *a) (1 + k2*a)*exp(-1im *k2 *a) ; (1im*k2)*exp(1im *k2 *a) (-1im*k2)*exp(-1im *k2 *a)]
    elseif k2 == 0
        print("k2 is zero")
        return [(exp(-1im * k1 * a)/2) ((k1*a - 1im)*exp(-1im * k1 * a)/(2k1)) ; (exp(1im * k1 * a)/2) ((k1*a + 1im)*exp(1im * k1 * a)/(2k1))]    
    else
        return (1/(2*k1))*[(k1+k2)*exp(-1im*a*(k1-k2)) (k1-k2)*exp(-1im*a*(k1+k2)); (k1-k2)*exp(1im*a*(k1+k2))  (k1+k2)*exp(1im*a*(k1-k2))]
    end
end

function getWavefunction(A, B, k, x)
    return A*exp.(1im.*x.*k).+B*exp.(-1im.*x.*k)
end

function getTotalWaveVector(E, U)
    return [getWaveVector(E, i) for i in U]
end

function getTransmissionProbability(t11)
    return (1/modulus(t11))^2
end

function getReflectionProbability(t21, t11)
    return (modulus(t21)^2)/(modulus(t11)^2)
end

""" N Reigon Simulations """

function nReigonTm(E, U, boundaries)
   
    k = getTotalWaveVector(E, U)
    TM = [1 0; 0 1]
    
    for i in 1:(length(boundaries))
        TM = TM * getTransferMatrix(k[i], k[i+1], boundaries[i])    
    end
    return TM
end


#=
function nReigonTm(E, U, boundaries)
    return getTransferMatrix(getWaveVector(E, U[1]), getWaveVector(E, U[2]), boundaries[1])   
end
=#

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

""" Energy Loops """
function energyLoop(E, U, boundaries)
    TM = [nReigonTm(i, U, boundaries) for i in E]
    
    t11, tp, rp = complex(zeros(length(E))), complex(zeros(length(E))), complex(zeros(length(E)))
    j = 1
    for i in TM
        t11[j] = i[1,1]
        tp[j] = getTransmissionProbability(i[1,1])
        rp[j] = getReflectionProbability(i[2,1], i[1,1])
        j = j + 1
    end
    
    #p1 = plot(real(tp))
    #p1 = plot!(real(rp))
    #display(plot(real(t11)))
    #display(p1)

    return t11
end

