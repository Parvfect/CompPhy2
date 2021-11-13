
include("potentials.jl")


parameters = [
    Bcn = [1.0, 0.0],
    U, An = piecewisePotential(0, 5, 11, gaussianPotential)
]

function gaussianFunction()
    # Scale x = xe-10
    new_gaussian = createGaussianPotential(5, 2.5, 1)
    Bcn = [1.0, 0.0]
    U, An = piecewisePotential(0, 5, 11, new_gaussian)

    transmissionArr = complex(zeros(0))
    reflectionArr = complex(zeros(0))
    E = 0:0.01:10
        
    @time for i in E
        grid, psi, tp, rp = nRegion(U, i, An, Bcn, 1e-2, 5)
        append!(transmissionArr, tp)
        append!(reflectionArr, rp)
    end

    plot(E, real(transmissionArr), title = "Transmission Probability", xlabel = "E", ylabel = "T(E)", label = "TP")
    #plot(E, real(reflectionArr), label = "RP")
end
