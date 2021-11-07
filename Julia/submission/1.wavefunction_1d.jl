# Generalizing schrodinger.jl for more robust simulations

using Plots
using LinearAlgebra

function Schrodinger1D(U, E, Bc2, a)
    #= U - Potentials of both states, E - Energy of Wavefunction, A2 & B2 - Boundary Conditions, m : mass
        Determines A1, B1 - First Boundary Reigon Coefficents given the paramters
        Uses Transfer Matrix Method =#
    
    # Initializing First Boundary Conditions
    A1 = 0 
    B1 = 0

    # Calculating wave vectors
    k1 = getWaveVector(E, U[1])
    k2 = getWaveVector(E, U[2])
    
    # Preparing matrices for transfer matrix method
    T = formTransferMatrix(k1, k2, a)
    Bc1 = [A1, B1]

    # Performing Matrix Multiplication and extracting values
    Bc1 = T * Bc2 

    return Bc1, k1, k2, T
end


function Ψ(A, B, k, x)
    return A*exp(im * k* x) + B * exp(- im * k* x) 
end

function probabilityAmplitude(v)
    #= Finds the Probability Amplitude of the Observable/Vector =#
    return dot(v, conj(v))
end
    
function getWaveVector(E, U)
    #= Returns the Wave Vector of the Electron given the Energy and the Potential of the State =#
    return sqrt(2* me *(Complex(E - U)))/hbar
end


function formTransferMatrix(k1, k2, a)
    #= Forming the Transfer Matrix for two generic boundaries given the wave vectors and the boundary location in the Grid=#

    k1k2 = k1 + k2
    k1_k2 = k1 - k2
    return 1/2k1 * [k1k2*exp(-1im*a*k1_k2) k1_k2*exp(-1im*a*k1k2) ; k1_k2*exp(1im*a*k1k2) k1k2*exp(1im*a*k1_k2)]
end
     
    
function transmissionProbability(t, k0, kn)
    #= Calculates the Transmission Probability given the transfer matrix, arbitary Potential Values, Energies and the Boundary Conditions =#
    if isinf(t) || k0 == 0 || isinf(kn)
        print("Breakout detected")
    end 
    a = modComplex(1/t) * modComplex(1/t) * (real(kn) /k0)
    return a
end

function modComplex(cn)
    #= Returns the Modulus of a Complex Number =#
    return sqrt(cn * conj(cn))
end


function reflectionProbability(t1, t2)
    #= Calculates the Reflection Probability given the transfer matrix, arbitary Potential Values, Energies and the Boundary Conditions =#
    return (abs(t2/t1)) * (abs(t2/t1))
end


function piecewisePotential(lowerLimit, upperLimit, n_reigons, potentialFunction)
    #= Returns a potential array of an n_reigon approximation to the provided potential function that can be fed into an electron simulation=#

    lengthReigon = (abs(lowerLimit) + upperLimit) / n_reigons
    U = zeros(n_reigons)
    ptr = lengthReigon/2

    for i in 1:n_reigons
        # Append the Potential Value for the midpoint of the ith reigon into the Potential Array
        U[i] =  potentialFunction(ptr)
        ptr += lengthReigon
    end

    return U, [lengthReigon for i in 1:n_reigons]
end

function generalisedWavefunction(grid, Bc, k)
    #= Calculates the Generalised Wavefunction for given reigon, Boundary Condition and Wave Vectors =#
    
    psi = complex(zeros(length(grid)))
    
    for i in 1:length(grid)
        psi[i] = Ψ(Bc[1], Bc[2], k, grid[i])
    end 

    return psi
end

function formGrid(start, step, limit)
    #= Forming the Grid for the Schrodinger Equation =#
    return start:step:limit
end

function nRegion(U, E, An, Bcn, step, lowerLimit = 0, upperLimit = 0)
    
    grid = zeros(0)
    psi = complex(zeros(0))

    if lowerLimit == 0 
        upperLimit = sum(An)
    end

    Bc2 = Bcn
    k1 = 0
    k2 = 0
    kn = 0
    t = 0
    tp = 0
    rp = 0

    for i in reverse(1:length(An))
        
        # Handling first grid of the System
        if i == 1
            grid_temp = formGrid(lowerLimit, step, upperLimit)
            psi_temp = generalisedWavefunction(grid_temp, Bc2, k1)
            prepend!(grid, grid_temp)
            prepend!(psi, psi_temp)
            tp =  transmissionProbability(t[1,1], k1, kn)
            rp = reflectionProbability(t[1,1], t[2,1])
            break
        end

        # For each step into the potential array we calculate reflection between two reigons
        boundary = upperLimit - An[i]
        grid_temp = formGrid(boundary, step, upperLimit)
        
        # Getting Transfer Matrix and Boundary Conditions for n and n-1th reigon 
        Bc1, k1, k2, t = Schrodinger1D(U[i-1 : i], E, Bcn, (upperLimit - An[i]))

        # Getting Wavefunction for the n-1th - n reigon 
        psi_temp = generalisedWavefunction(grid_temp, Bc2, k2)
        
        # Append the values of the reigon to the start of the main grid
        prepend!(grid, grid_temp)
        prepend!(psi, psi_temp)

        upperLimit = boundary
        Bc2 = Bc1


        if i == length(An)
            kn = k2
        end
    end

    return grid, psi, tp, rp
end

function plotSimulation(grid, psi, energy)
    #= Plots the Simulation of the Wavefunction through N Reigon Potential =#
    
    display(plot(grid, real(psi), title = "Wavefunction in one dimension, E = $energy", label = "Real part"))
    xlabel!("Positon in grid (x)")
    ylabel!("Wavefunction")
    plot!(grid, imag(psi), label = "Imaginary part")
end


