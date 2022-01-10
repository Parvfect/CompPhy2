#= ../harmonic_oscillator.jl
System Parameters for the Harmonic Oscillator 
=#

e, me, hbar, A3 = 1.6e-19, 9.11e-31, 1.05e-34, 1.0
save_path = "C:/Users/Parv/Documents/compphy/Julia/Data/"

function createHarmonicPotential(w)
    return x -> (x^2*w^2*me)/2
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

p = createHarmonicPotential(1.5e14)
U, an = harmonicPiecewise(10, 11, p)
U, reigon_lengths = U*e, an*1e-9
boundaries = getBoundaries(an*1e-9)
