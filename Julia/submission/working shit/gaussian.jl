using Plots
using LinearAlgebra
include("generic_soln.jl")

e, me, hbar, A3 = 1.6e-19, 9.11e-31, 1.05e-34, 1.0
save_path = "C:/Users/Parv/Documents/compphy/Julia/Data/"
#U, reigon_lengths, boundaries = [3, 0, 3]*e, [2e-9, 8e-9, 2e-9], [2e-9, 10e-9]


function createGaussianPotential(V, a)
    #= V = 3 , a = 5 =#
    return x -> V * exp(- (x-a)^2/6.25)
end


function piecewisePotential(upperLimit, n_reigons, potentialFunction)
    lengthReigon = upperLimit / n_reigons
    U = zeros(n_reigons)
    ptr = lengthReigon/2

    for i in 1:n_reigons
        U[i] =  potentialFunction(ptr)
        ptr += lengthReigon
    end

    return U, [lengthReigon for i in 1:n_reigons]
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
energyLoop(E, U, boundaries)

#tprpSim()
#nReigonSim(0.75*e, U, boundaries, an*1e-9)
# Bound States - Highest to Lowest
#E = 3.843359e-19
#E = 1.2678747928e-19
#E = 7.1410000239e-20
#E = 1.4122357142399997e-20
#E = 3.531310235e-21
#E = 0.8828720751e-21

