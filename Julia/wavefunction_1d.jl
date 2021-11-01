# Generalizing schrodinger.jl for more robust simulations

using Plots
using LinearAlgebra

# Constants
e=1.6e-19
me=9.11e-31
# a=3e-10
hbar=1.05e-34



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

    return Bc1, k1, k2
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
     
    
function transmissionProbability(t, k1, kn)
    #= Calculates the Transmission Probability given the transfer matrix, arbitary Potential Values, Energies and the Boundary Conditions =#
    return (modulus(1/ t[1,1])) * (modulus(1/ t[1,1])) * (kn / k1)

end

function reflectionProbability(t, U, E, A2, B2)
    #= Calculates the Reflection Probability given the transfer matrix, arbitary Potential Values, Energies and the Boundary Conditions =#
    return (modulus(t[2,1]/t[1,1])) * (modulus(t[2,1]/t[1,1]))
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

function nRegion(U, E, An, Bcn, step)
    
    grid = zeros(0)
    psi = complex(zeros(0))

    upperLimit = sum(An)
    Bc2 = Bcn
    k1 = 0
    k2 = 0

    for i in reverse(1:length(An))

        if i == 1
            grid_temp = formGrid(0, step, upperLimit)
            psi_temp = generalisedWavefunction(grid_temp, Bc2, k1)
            prepend!(grid, grid_temp)
            prepend!(psi, psi_temp)
            break
        end

        # For each step into the potential array we calculate reflection between two reigons
        boundary = upperLimit - An[i]
        #lowerLimit = boundary - An[i-1]

        grid_temp = formGrid(boundary, step, upperLimit)
        #grid_1 = formGrid(lowerLimit, step, boundary)

        # Getting Transfer Matrix and Boundary Conditions for n and n-1th reigon 
        Bc1, k1, k2 = Schrodinger1D(U[i-1 : i], E, Bcn, (upperLimit - An[i]))

        # Getting Wavefunction for the n-1th - n reigon 
        psi_temp = generalisedWavefunction(grid_temp, Bc2, k2)
        #psi_1 = generalisedWavefunction(grid_1, Bc1, k1)
        
        # Found a problem with the waefunction 
        # the wavefcuntions in the middle won't be appended directly only the first and last would, the 
        # others would just have a changed wave vector and boundary conditions

        # Append the values of the reigon to the start of the main grid
        prepend!(grid, grid_temp)
        prepend!(psi, psi_temp)

        upperLimit = boundary
        Bc2 = Bc1

    end

    return grid, psi
end

function plotSimulation(grid, psi)
    #= Plots the Simulation of the Wavefunction through N Reigon Potential =#
    
    display(plot(grid, real(psi)))
    plot!(grid, imag(psi))
end


function default_simulation()
    U = [2,0,2] * e
    E = [0.75, 1.5, 2.5] * e
    Bc2 = [1.0, 0.0]
    size_reigon = 3e-10 + 2e-9
    An = [size_reigon, size_reigon, size_reigon]

    grid, psi = nRegion(U, E[1], An, Bc2, 1e-11)
    plotSimulation(grid, psi)
end

default_simulation()