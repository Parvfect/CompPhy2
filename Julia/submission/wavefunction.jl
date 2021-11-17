using Plots
using LinearAlgebra

include("constants.jl")

""" Helper Methods """

function Ψ(A, B, k, x)
    return A*exp(1im * k* x) + B * exp(- 1im * k* x) 
end

function probabilityAmplitude(v)
    return dot(v, conj(v))
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

function modComplex(cn)
    return sqrt(cn * conj(cn))
end

function reflectionProbability(t1, t2)
    return (abs(t2/t1)) * (abs(t2/t1))
end

function formGrid(start, step, limit)
    return start:step:limit
end


""" Integration Methods """

function transferMatrixMethod(U, E, Bc2, a)
   
    A1 = 0 
    B1 = 0

    k1 = getWaveVector(E, U[1])
    k2 = getWaveVector(E, U[2])
    T = formTransferMatrix(k1, k2, a)
    Bc1 = [A1, B1]

    Bc1 = T * Bc2 

    return Bc1, k1, k2, T
end

function generalisedWavefunction(grid, Bc, k)
    psi = complex(zeros(length(grid)))
    
    for i in 1:length(grid)
        psi[i] = Ψ(Bc[1], Bc[2], k, grid[i])
    end 

    return psi
end


""" Plotting functions """

function plotWavefunction(grid, psi, energy)
    
    display(plot(grid, real(psi), title = "Wavefunction in one dimension, E = $energy", label = "Real part"))
    xlabel!("Positon in grid (x)")
    ylabel!("Wavefunction")
    plot!(grid, imag(psi), label = "Imaginary part")
end

function plotTprp(tp, rp, energy)
        
    display(plot(energy, tp, label = "Transmission probability", xlabel = "Energy", ylabel = "Probability"))
    xlabel!("Energy")
    ylabel!("Probability")
    plot!(energy, rp, label = "Reflection probability")
end


""" n reigon simulations """

function t11Sim(U, E, An, Bcn, step)
    
    grid = zeros(0)
    psi = complex(zeros(0))
    Bc2, Bc1, k1, k2, t, kn = Bcn, 0, 0, 0, 0, 0

    upperLimit = sum(An)
    
    for i in reverse(1:length(An))
        
        boundary = upperLimit - An[i]
        grid_temp = formGrid(boundary, step, upperLimit)

        if i != 1
            Bc1, k1, k2, t = transferMatrixMethod(U[i-1 : i], E, Bc2, boundary)
        else
            #println(t[1,1], k1, kn)
            return t[1,1], transmissionProbability(t[1,1], k1, kn), reflectionProbability(t[1,1], t[2,2])
        end
        
        if i== length(An)
            kn = k1
        end

        upperLimit = boundary
        Bc2 = Bc1
    end

end

function psiSim(U, E, An, Bcn, step)
    
    grid = zeros(0)
    psi = complex(zeros(0))
    Bc2, Bc1, k1, k2 = Bcn, 0, 0, 0, 0

    upperLimit = sum(An)

    for i in reverse(1:length(An))
        
        boundary = upperLimit - An[i]
        grid_temp = formGrid(boundary, step, upperLimit)

        if i != 1
            Bc1, k1, k2, t = transferMatrixMethod(U[i-1 : i], E, Bc2, (upperLimit - An[i]))
            psi_temp = generalisedWavefunction(grid_temp, Bc2, k2)  
        else
            psi_temp = generalisedWavefunction(grid_temp, Bc2, k1)
        end

        if i == 100
            println("Suck a dic bitch ", k1)
        end
        prepend!(grid, grid_temp)
        prepend!(psi, psi_temp)

        upperLimit = boundary
        Bc2 = Bc1
    end
    
    return grid, psi
end

