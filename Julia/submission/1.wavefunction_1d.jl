using Plots
using LinearAlgebra

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


function piecewisePotential(lowerLimit, upperLimit, n_reigons, potentialFunction)
    lengthReigon = (abs(lowerLimit) + upperLimit) / n_reigons
    U = zeros(n_reigons)
    ptr = lengthReigon/2

    for i in 1:n_reigons
        U[i] =  potentialFunction(ptr)
        ptr += lengthReigon
    end

    return U, [lengthReigon for i in 1:n_reigons]
end

function generalisedWavefunction(grid, Bc, k)
    psi = complex(zeros(length(grid)))
    
    for i in 1:length(grid)
        psi[i] = Ψ(Bc[1], Bc[2], k, grid[i])
    end 

    return psi
end


function nRegion(U, E, An, Bcn, step, upperLimit = 0)
    
    grid = zeros(0)
    psi = complex(zeros(0))
    Bc2, Bc1, k1, k2, kn, t, tp, rp = Bcn, 0, 0, 0, 0, 0, 0, 0

    if upperLimit == 0 
        upperLimit = sum(An)
    end

    for i in reverse(1:length(An))
        
        boundary = upperLimit - An[i]
        grid_temp = formGrid(boundary, step, upperLimit)

        if i != 1
            Bc1, k1, k2, t = transferMatrixMethod(U[i-1 : i], E, Bc2, (upperLimit - An[i]))
            psi_temp = generalisedWavefunction(grid_temp, Bc2, k2)  
        else
            psi_temp = generalisedWavefunction(grid_temp, Bc2, k1)
            tp =  transmissionProbability(t[1,1], k1, kn)
            rp = reflectionProbability(t[1,1], t[2,1])    
        end
        
        if i == length(An)
            kn = k2
        end

        prepend!(grid, grid_temp)
        prepend!(psi, psi_temp)

        upperLimit = boundary
        Bc2 = Bc1
    end
    
    return grid, psi, tp, rp
end

function plotSimulation(grid, psi, energy)
    
    display(plot(grid, real(psi), title = "Wavefunction in one dimension, E = $energy", label = "Real part"))
    xlabel!("Positon in grid (x)")
    ylabel!("Wavefunction")
    plot!(grid, imag(psi), label = "Imaginary part")
end


""" Helper Methods """


function Ψ(A, B, k, x)
    return A*exp(im * k* x) + B * exp(- im * k* x) 
end

function probabilityAmplitude(v)
    return dot(v, conj(v))
end
    
function getWaveVector(E, U)
    return sqrt(2* me *(Complex(E - U)))/hbar
end


function formTransferMatrix(k1, k2, a)
    k1k2 = k1 + k2
    k1_k2 = k1 - k2
    return 1/2k1 * [k1k2*exp(-1im*a*k1_k2) k1_k2*exp(-1im*a*k1k2) ; k1_k2*exp(1im*a*k1k2) k1k2*exp(1im*a*k1_k2)]
end
     
    
function transmissionProbability(t, k0, kn)
    if isinf(t) || k0 == 0 || isinf(kn)
        print("Breakout detected")
    end 
    a = modComplex(1/t) * modComplex(1/t) * (real(kn) /k0)
    return a
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
