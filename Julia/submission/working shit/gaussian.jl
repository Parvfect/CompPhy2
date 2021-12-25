include("generic_soln.jl")

# Transmission Probability Plot looks sound - 25/12/21

e, me, hbar, A3 = 1.6e-19, 9.11e-31, 1.05e-34, 1.0
save_path = "C:/Users/Parv/Documents/compphy/Julia/Data/"
#U, reigon_lengths, boundaries = [3, 0, 3]*e, [2e-9, 8e-9, 2e-9], [2e-9, 10e-9]


function createGaussianPotential(V, a)
    #= V = 3 , a = 5 =#
    return x -> V * exp(- (x-a)^2/6.25)
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

p = createGaussianPotential(2,5)
U, an = piecewisePotential(10, 1111, p)
U = U*e
boundaries = getBoundaries(an*1e-9)
E = 1e-21:1e-22:10e-19
t11 = energyLoop(E, U, boundaries)
#display(plot(real(t11)))
