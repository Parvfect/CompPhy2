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
    
    plotWaveFunction1D(Bc1, Bc2, k1, k2)

    return T, Bc1, k1, k2
end


function Ψ(A, B, k, x)
    return A*exp(im * k* x) + B * exp(- im * k* x) 
end

function probabilityAmplitude(v)
    #= Finds the Probability Amplitude of the Observable/Vector =#
    return dot(v, conj(v))
end

function plotWaveFunction1D(Bc1, Bc2, k1, k2)
    # Plots the Wavefunction of an electron incident on the boundary in 1d """ 
    dx = 0.02
    grid = -2e-9:1e-11:2e-9
    psi = complex(zeros(length(grid)))
    
    for i in 1:length(grid)

        if grid[i] < a
            psi[i] = Ψ(Bc1[1], Bc1[2], k1, grid[i])
        
        # Boundary in the middle 
        else
            psi[i] = Ψ(Bc2[1], Bc2[2], k2, grid[i])
        end
    end

    plot(grid, real(psi),  label = "Wavefunction for Time Independent Schordinger Equation in One Dimension", xlabel = "x", ylabel = "psi(x)")
    plot!(grid, imag(psi))
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

function generalisedWavefunction(grid, Bc1, Bc2, k1, k2)
    #= Calculates the Generalised Wavefunction for given reigon, Boundary Condition and Wave Vectors =#
    
    psi = complex(zeros(length(grid)))
    
    for i in 1:length(grid)

        if grid[i] < a
            psi[i] = Ψ(Bc1[1], Bc1[2], k1, grid[i])
        
        # Boundary in the middle 
        else
            psi[i] = Ψ(Bc2[1], Bc2[2], k2, grid[i])
        end
    end

    return psi
end

function formGrid(start, limit, step)
    #= Forming the Grid for the Schrodinger Equation =#
    return start:step:limit
end

function nRegion(U, E, An, Bcn, lowerLimit, upperLimit, step)
    #= Creates and Combines n reigons for Wavefunction Simulation
        U - Array of Potentials of a each reigon
        E - Energy of State
        An - Array of Boundary locations
        Bcn - Array of Boundary Conditions for the nth reigon 

    Working   
    1. Gets Transfer Matrix and Boundary Conditions for n-1th - nth reigon 
    2. Gets Wavefunction for the n-1th - n reigon
    3. Updates the reigon boundaries and the boundary conditions for the next iteration  
    4. Repeat 1-3 until U is empty
    5. Return the grid and wavefunctions array

    Assumes all reigons are symmetric - needs to be changed for general case
    =#
    
    grid = formGrid(lowerLimit, upperLimit, step)
    lowerLimit = An[length(An)] - upperLimit
    Wavefunctions = Complex(zeros(length(grid)))

    for i in length(U):-1:2
        
        # Getting Transfer Matrix and Boundary Conditions for n-1th reigon 
        T, Bc, k1, k2 = Schrodinger1D(U[i, i - 1], E, Bcn, A[i])

        # Getting Wavefunction for the n-1th - n reigon 
        psi[lowerLimit: upperLimit] = generalisedWavefunction(grid[lowerLimit: upperLimit], Bc1, Bc2, k1, k2)
        
        # Updating Boundary Condition
        Bcn = Bc

        # Defining the boundaries of the n-1th reigon - symmetric assumption
        upperLimit = lowerLimit
        lowerLimit = An[i - 1] - (upperLimit - An[i])

    end

    return grid, Wavefunctions
end

function plotSimulation()
    #= Plots the Simulation of the Wavefunction through N Reigon Potential =#
end

U = [0, 2] * e
E = [0.75, 1.5, 2.5] * e
Bc2 = [1.0, 0.0]

T, Bc1, k1, k2 = Schrodinger1D(U, E[1], Bc2, 3e-10)
plotWaveFunction1D(Bc1, Bc2, k1, k2)