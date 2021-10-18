using Plots
using LinearAlgebra


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
    T = 1/2k1 * 
    [k1k2*exp(-im*k1_k2) k1_k2*exp(-im*k1k2) ; k1_k2*exp(im*k1k2) k1k2*exp(im*k1_k2)]
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
    grid = 0:dx:10
    psi = complex(zeros(length(grid)))
    
    for i in 1:length(grid)

        if grid[i] < 5
            psi[i] = Ψ(Bc1[1], Bc1[2], k1, grid[i])
        
        # Boundary in the middle 
        else
            psi[i] = Ψ(Bc2[1], Bc2[2], k2, grid[i])
        end
    end

    mod_psi = [probabilityAmplitude(i) for i in psi]
    plot(grid, real(psi),  label = "Wavefunction for Time Independent Schordinger Equation in One Dimension", xlabel = "x", ylabel = "psi(x)")
    #plot!(grid, real(psi))
    #plot!()
end
    
    

U = [0, 2]
E = [0.5, 1.5, 2.5]
Schrodinger1D(U, E[1], 2.4, 0)