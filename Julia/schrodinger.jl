using Plots
using LinearAlgebra

e = MathConstants.e

function Schrodinger1D(U, E, A2, B2)
    #= U - Potentials of both states, E - Energy of Wavefunction, A2 & B2 - Boundary Conditions, m : mass
        Determines A1, B1 - First Boundary Reigon Coefficents given the paramters
        Uses Transfer Matrix Method =#
    
    # Initializing First Boundary Conditions
    A1 = 0 
    B1 = 0

    # Calculating constants
    k1 = sqrt(complex(2* (E - U[1])))/6.626
    k2 = sqrt(complex(2* (E - U[2])))/6.626
    k1k2 = k1 + k2
    k1_k2 = k1 - k2

    # Forming matrix
    T = [k1k2*e^(-im*k1_k2) k1_k2*e^(-im*k1k2) ; k1_k2*e^(im*k1k2) k1k2*e^(im*k1_k2)]
    Bc2 = [A2, B2] / 2k1
    Bc1 = [A1, B1]

    # Performing Matrix Multiplication and extracting values
    Bc1 = T * Bc2  
    A1 = Bc1[1]
    B1 = Bc1[2]

    plotWaveFunction1D(Bc1, Bc2, k1, k2)

end

#=
function Ψ()
end
=#

function plotWaveFunction1D(Bc1, Bc2, k1, k2)
    # Plots the Wavefunction of an electron incident on the boundary in 1d """ 
    dx = 0.02
    grid1 = 0:dx:5
    grid2 = 5:dx:10
    psi1 = []
    psi2 = []
    a = 100.0
    
    append!(psi1, Bc1[1]e^(im*k1*grid1) + Bc1[2]e^(-im*k1*i))
    append!(psi2, Bc2[1]e^(im*k2*grid2) + Bc2[2]e^(-im*k2*i))

    plot(Complex(grid1), psi1, label = "Wavefunction for Time Independent Schordinger Equation in One Dimension", xlabel = "x", ylabel = "psi(x)")
    plot!(Complex(grid2), psi2)
end
    
    

U = [0, 2]
E = [0.5, 1.5, 2.5]
Schrodinger1D(U, E[1], 2.4, 0)