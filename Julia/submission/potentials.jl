
w = 1.5e14
m = 9.1e-31

using Statistics

function createGaussianPotential(V, a)
    #= V = 3 , a = 5 =#
    return x -> V * exp(- (x-a)^2/6.25)
end

        
function createLinePotential(m, c) 
    return x -> (m*x + c)
end

function createHarmonicPotential()
    return x -> (x^2*w^2*m)/2
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

p = createGaussianPotential(3,5)
U, an = piecewisePotential(10, 5, p)
