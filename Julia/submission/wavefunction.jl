using Plots
using LinearAlgebra

include("constants.jl")

""" Helper Methods """

function Ψ(A, B, k, x)
    return A*exp(1im * k* x) + B * exp(- 1im * k* x) 
end

function getWaveVector(E, U)
    return (sqrt(2* me *(Complex(E - U)))/hbar)e-9
end

function formTransferMatrix(k1, k2, a)
    if k1 == 0
        return [(1 - k2*x)*exp(1im *k2 *x) (1 + k2*x)*exp(-1im *k2 *x) ; (1im*k2)*exp(1im *k2 *x) (-1im*k2)*exp(-1im *k2 *x)]
    elseif k2 == 0
        return [(exp(-1im * k1 * x)/2) ((k1*x - 1im)*exp(-1im * k1 * x)/(2k*1)) ; (exp(1im * k1 * x)/2) ((k1*x + 1im)*exp(1im * k1 * x)/(2k*1))]    
    else
        return  1/2k1 * [(k1 + k2)*exp(-1im*a*(k1 - k2)) (k1 - k2)*exp(-1im*a*(k1 + k2)) ; (k1 - k2)*exp(1im*a*(k1 + k2)) (k1 + k2)*exp(1im*a*(k1 - k2))]
    end
end 
    
function transmissionProbability(t, k0, kn)
    return abs(1/t) * abs(1/t) * ((kn) /k0)
end

function reflectionProbability(t1, t2)
    return (abs(t2/t1)) * (abs(t2/t1))
end

""" Integration Methods """

function transferMatrixMethod(k1, k2, Bc2, a)
    T = formTransferMatrix(k1, k2, a)
    return T, T * Bc2 
end

function totalWavefunction(An, A, B, K)
    psi, ptr = [], 0
    for i in 1:length(An)
        an, a, b, k = An[i], A[i], B[i], K[i]
        grid_temp = ptr:0.01:(ptr + an)
        psi = append!(psi, generalisedWavefunction(a, b, k, grid_temp))
        ptr += an
    end
    display(plot(real(psi)))
    return psi 
end

function generalisedWavefunction(A, B, K, grid)
    return [Ψ(A, B, K, i) for i in grid]    
end


""" Plotting functions """

function plotWavefunction(grid, psi, energy)
    display(plot(grid, real(psi), title = "Wavefunction in one dimension, E = $energy", label = "Real part"))
    xlabel!("Positon in grid (x)")
    ylabel!("Wavefunction")
    #plot!(grid, imag(psi), label = "Imaginary part")
end

function plotTprp(tp, rp, energy)
    display(plot(energy, tp, label = "Transmission probability", xlabel = "Energy", ylabel = "Probability"))
    xlabel!("Energy")
    ylabel!("Probability")
    plot!(energy, rp, label = "Reflection probability")
end

function plotT11(t11, energy)
    display(plot(energy, t11, label = "T11", xlabel = "Energy", ylabel = "Variation of T[1,1] over energy range"))
    xlabel!("Energy")
    ylabel!("T[1,1]")
end

""" n reigon simulations """

function solveTMM(U, E, An, Bcn, step)
    
    k = [getWaveVector(E, i) for i in U]
    A, B = complex(zeros(length(U))), complex(zeros(length(U)))
    T, Bc2 = 0, Bcn
    upperLimit = sum(An)
    
    for i in reverse(1:length(An))
        A[i], B[i] = Bc2[1], Bc2[2]
        if i == 1
            return k, A, B, T[1,1] # Would need to be changed for rp
        end
        T, Bc1 = transferMatrixMethod(k[i], k[i-1], Bc2, upperLimit-An[i])
        upperLimit -= An[i]
        Bc2 = Bc1
    end
end

function energyLoop(E, U, An, step)
    
    t11arr, k_arr, A_arr, B_arr = zeros(0), complex(zeros(0)), complex(zeros(0)), complex(zeros(0))
    for i in E
        K, A, B, t11 = solveTMM(U, i, An, [1.0, 0.0], 1e-2)
        t11arr = append!(t11arr, real(t11))
        k_arr = append!(k_arr, K)
        A_arr = append!(A_arr, A)
        B_arr = append!(B_arr, B)
    end
    return t11arr, k_arr, A_arr, B_arr
end
