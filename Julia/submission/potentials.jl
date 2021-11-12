
w = 1.5e14
m = 9.1e-31


function createGaussianPotential(a, b, c)
    return x -> a * exp( - ((x - b)^2) / (2*c^2))
end

function createLinePotential(m, c) 
    return x -> (m*x + c)
end

function createHarmonicPotential()
    return x -> (x^2*w^2*m)/2
end