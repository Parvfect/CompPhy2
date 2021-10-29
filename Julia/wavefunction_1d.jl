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

    return T, Bc1, k1, k2
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
    return 1/2k1 * 
    [(k1 + k2)*exp(-1im*a*(k1 - k2)) (k1 - k2)*exp(-1im*a*(k1 + k2)) ; (k1 - k2)*exp(1im*a*(k1 + k2)) (k1 + k2)*exp(1im*a*(k1 - k2))]  
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

function nRegion(U, E, An, Bcn, lowerLimit, upperLimit, step)
    #= Creates and Combines n reigons for Wavefunction Simulation
        U - Array of Potentials of a each reigon
        E - Energy of State
        An - Array of lengths of reigons
        Bcn - Array of Boundary Conditions for the nth reigon 

    Working   
    1. Gets Transfer Matrix and Boundary Conditions for n-1th - nth reigon 
    2. Gets Wavefunction for the n-1th - n reigon
    3. Updates the reigon boundaries and the boundary conditions for the next iteration  
    4. Repeat 1-3 until U is empty
    5. Return the grid and wavefunctions array

    Assumes all reigons are symmetric - needs to be changed for general case
    =#
    
    grid = zeros(0)
    psi = complex(zeros(0))

    upperLimit = sum(An)
    Bc2 = Bcn

    for i in reverse(1:length(U))

        if i == 1
            break
        end
        # For each step into the potential array we calculate reflection between two reigons
        boundary = upperLimit - An[i]
        lowerLimit = boundary - An[i-1]

        grid_2 = formGrid(upperLimit, (-step), boundary)
        print(grid_2)
        grid_1 = formGrid(boundary, (-step), lowerLimit)

        # Getting Transfer Matrix and Boundary Conditions for n and n-1th reigon 
        T, Bc1, k1, k2 = Schrodinger1D(U[i-1 : i], E, Bcn, (upperLimit - An[i]))

        # Getting Wavefunction for the n-1th - n reigon 
        psi_2 = generalisedWavefunction(grid_2, Bc2, k1)
        psi_1 = generalisedWavefunction(grid_1, Bc1, k2)
        
        # Append the values of the reigon to the start of the main grid
        prepend!(grid, grid_2, grid_1)
        prepend!(psi, psi_2, psi_1)

        # Defining the boundaries of the n-1th reigon - symmetric assumption
        grid_temp = []
        Bc2 = Bc1

        upperLimit = lowerLimit

    end

    return grid, psi
end

function plotSimulation(grid, psi)
    #= Plots the Simulation of the Wavefunction through N Reigon Potential =#
    print(typeof(real(psi)))
    #display(plot(grid, real(psi)))
    #plot!(grid, imag(psi))
end


function default_simulation()
    U = [0, 2] * e
    E = [0.75, 1.5, 2.5] * e
    Bc2 = [1.0, 0.0]
    a = 1e-10

    T, Bc1, k1, k2 = Schrodinger1D(U, E[1], Bc2, a)
    plotWaveFunction1D(Bc1, Bc2, k1, k2, a)
end


default_simulation()