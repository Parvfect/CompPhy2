#=
25/12/21
Going out of memory for the bound states
Sign Flip loss for the other bound states
t11 looks solid
=#


using Plots
using LinearAlgebra
include("generic_soln.jl")
include("bound_states.jl")

e, me, hbar, A3 = 1.6e-19, 9.11e-31, 1.05e-34, 1.0
save_path = "C:/Users/Parv/Documents/compphy/Julia/Data/"
#U, reigon_lengths, boundaries = [3, 0, 3]*e, [2e-9, 8e-9, 2e-9], [2e-9, 10e-9]

function createHarmonicPotential(w)
    return x -> (x^2*w^2*me)/2
end

function nReigonSim(E, U, boundaries, reigon_lengths)
    A, B = nReigon(E, U, 1.0, 0, boundaries)
    nReigonPlot(A, B, getTotalWaveVector(E, U), reigon_lengths, false)
end


function t11Sim(U, boundaries)
    E = 1e-22:1e-22:3e-19
    t11 = energyLoop(E, U, boundaries)
    display(plot(real(t11)))
end

function getBoundaries(reigon_lengths)
    ptr = 0
    boundaries = zeros(length(reigon_lengths) - 1)
    for i in 1:(length(reigon_lengths)-1)
        boundaries[i] = ptr + reigon_lengths[i]
        ptr = ptr + reigon_lengths[i]
    end
    return boundaries
end

function plotBoundStates(U, boundaries)
    # Still got repeats in bound states, and its not picking up all bound states in the reigon
    E = 1e-22:1e-22:3e-19
    boundStates = getAllBoundStates(E, U, boundaries)
end


p = createHarmonicPotential(1.5e14)
U, an = harmonicPiecewise(10, 11, p)
display(plot(U))
U = U*e
boundaries = getBoundaries(an*1e-9)
#t11Sim(U, boundaries)
plotBoundStates(U, boundaries)
#nReigonSim(9.557812500000001e-21, U, boundaries, an*1e-9)