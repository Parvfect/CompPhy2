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
    
    p1 = plot(real(tp))
    #p1 = plot!(real(rp))
    #display(plot(real(t11)))
    display(p1)

end

#function energyLoop(E, U, boundaries)
    
    #t11 = [nReigonTm(i, U, boundaries)[1,1] for i in E]
    
    #=t11, tp, rp = complex(zeros(length(E))), zeros(length(E)), zeros(length(E))
    j = 1
    for i in E
        TM = nReigonTm(i, U, boundaries)
        t11[j] = TM[1,1]
        #tp[j] = getTransmissionProbability(t11[i])
        #rp[j] = getReflectionProbability(TM[2,1], t11[i])
        j += 1
    end

    print(t11)
    =#
    #p1 = plot(t11, xlabel = "E", ylabel = "Transmission Probability", xlim=(-Inf,+Inf), ylim=(0,1), grid = true)
    #p1 = plot!(tp, xlabel = "E", ylabel = "Transmission Probability", xlim=(-Inf,+Inf), ylim=(0,1), grid = true)
    #p1 = plot!(rp, xlabel = "E", ylabel = "Reflection Probability", xlim=(-Inf,+Inf), ylim=(0,1), grid = true)
    #display(p1)
#end

e, me, hbar, A3 = 1.6e-19, 9.11e-31, 1.05e-34, 1.0
U = [0, 2, 0]*e
#E=1.4122357142399997e-20
reigon_lengths = [2e-9, 5e-9, 2e-9]
boundaries = [2e-9, 7e-9]
E = 1e-22:1e-22:6e-19 
t11 = energyLoop(E, U, boundaries)
#display(plot(real(t11)))