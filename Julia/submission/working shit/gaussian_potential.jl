#= ../gaussian_potential.jl
System Parameters for the Gaussian Potential 
=#

e, me, hbar, A3 = 1.6e-19, 9.11e-31, 1.05e-34, 1.0
save_path = "C:/Users/Parv/Documents/compphy/Julia/Data/gaussian_potential/"


function createGaussianPotential(V, a)
    #= V = 3 , a = 5 =#
    return x -> V * exp(-(x-a)^2)
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

p = createGaussianPotential(1,5)
U, an = piecewisePotential(10, 1111, p)
U, reigon_lengths = U*e, an*1e-9
boundaries = getBoundaries(reigon_lengths)
