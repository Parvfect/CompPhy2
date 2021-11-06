include("1.wavefunction_1d.jl")
include("potentials.jl")



# Constants
e=1.6e-19
me=9.11e-31
# a=3e-10
hbar=1.05e-34




function finiteWell()
    U = [3,0,3] * e
    E = [0.75, 1.5, 2.5] * e
    Bc2 = [1.0, 0.0]
    size_reigon = 3e-10 + 2e-9
    An = [size_reigon, U[1], size_reigon]

    grid, psi = @time nRegion(U, 20*e, An, Bc2, 1e-11)
    plotSimulation(grid, psi, 20*e)
end

function gaussianFunction()
    # Scale x = xe-10
    new_gaussian = createGaussianPotential(5*e*(10e9), 2.5, 1)
    Bcn = [1.0, 0.0]
    U, An = piecewisePotential(-2, 2, 1001, new_gaussian)
    #grid, psi, tp = @time nRegion(U, 0.02, An, Bcn, 1e-3)
    transmissionArr = complex(zeros(0))

    E = 1*e*(10e9):1e-11:20*e*(10e9)
    
    try
        @time for i in E
            grid, psi, tp = nRegion(U, i, An, Bcn, 1e-3)
            append!(transmissionArr, tp)
        end
    catch   
        print("Exception $i")
    end

    plot(E, real(transmissionArr), title = "Transmission Probability", xlabel = "E", ylabel = "T(E)")
end

gaussianFunction()

