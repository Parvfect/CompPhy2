using Plots
using LinearAlgebra

# Constants
e=1.6e-19
me=9.11e-31
a=3e-10
hbar=1.05e-34


function Schrodinger1D(U, E, A2, B2)
    #= U - Potentials of both states, E - Energy of Wavefunction, A2 & B2 - Boundary Conditions, m : mass
        Determines A1, B1 - First Boundary Reigon Coefficents given the paramters
        Uses Transfer Matrix Method =#
    
    # Initializing First Boundary Conditions
    A1 = 0 
    B1 = 0

    # Calculating constants
    k1 = sqrt(2* me *complex((E - U[1])))/hbar
    k2 = sqrt(2 * me *complex((E - U[2])))/hbar
   
    k1k2 = k1 + k2
    k1_k2 = k1 - k2

    # Forming matrix
    T = 1/2k1 * 
    [k1k2*exp(-1im*a*k1_k2) k1_k2*exp(-1im*a*k1k2) ; k1_k2*exp(1im*a*k1k2) k1k2*exp(1im*a*k1_k2)]
    Bc2 = [A2, B2]
    Bc1 = [A1, B1]

    # Performing Matrix Multiplication and extracting values
    Bc1 = T * Bc2 
    A1 = Bc1[1]
    B1 = Bc1[2]

    plotWaveFunction1D(Bc1, Bc2, k1, k2)

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
    
    

U = [0, 2] * e
E = [0.75, 1.5, 2.5] * e
Schrodinger1D(U, E[1], 1.0, 0)