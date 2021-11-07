

function createGaussianPotential(a, b, c)
    return x -> a * exp( - ((x - b)^2) / (2*c^2))
end

function createLinePotential(m, c) 
    return x -> (m*x + c)
end