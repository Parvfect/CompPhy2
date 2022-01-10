#= ../triangular_potential.jl
System Parameters for the Triangular Potential 
=#

save_path = "C:/Users/Parv/Documents/compphy/Julia/Data/triangular_potential/"
include("generic_soln.jl")

function triangularPotential(x)
    if x!=0 
        return 2 - abs(x)
    else
        return 0
    end
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
        

U, an = harmonicPiecewise(5, 111, triangularPotential)
U, reigon_lengths = U*1.6e-19, an*1e-10
boundaries = getBoundaries(reigon_lengths)
