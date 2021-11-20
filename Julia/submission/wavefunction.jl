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
    k1k2 = k1 + k2
    k1_k2 = k1 - k2
    t =  1/2k1 * [k1k2*exp(-1im*a*k1_k2) k1_k2*exp(-1im*a*k1k2) ; k1_k2*exp(1im*a*k1k2) k1k2*exp(1im*a*k1_k2)]
    return t
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

function generalisedWavefunction(grid, A, B, k)
    psi = complex(zeros(length(grid)))    
    for i in 1:length(grid)
        psi[i] = Ψ(A, B, k, grid[i])
    end 
    return psi
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
    A, B = zeros(length(U)), zeros(length(U))
    T, Bc2 = 0, Bcn
    upperLimit = sum(An)
    
    for i in reverse(1:length(An))
        A[i], B[i] = Bc2[1], Bc2[2]
        if i == 1
            return k, A, B, T[1,1] # Would need to be changed for rp
        end
        Bc1, k1, k2, t = transferMatrixMethod(k[i], k[i-1], Bc2, upperLimit-An[i])
        upperLimit -= An[i]
        Bc2 = Bc1
    end
end

function energyLoop(E, U, An, step)
    
    t11arr, k_arr, A_arr, B_arr = zeros(0), zeros(0), zeros(0), zeros(0)
    for i in E
        k, A, B, t11 = solveTMM(U, i, An, [1.0, 0.0], 1e-2)
        t11arr = append!(t11arr, real(t11))
        k_arr = append!(k_arr, K)
        A_arr = append!(A_arr, A)
        B_arr = append!(B_arr, B)
    end
    return t11arr, k_arr, A_arr, B_arr
end
