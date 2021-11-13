include("wavefunction.jl")
include("potentials.jl")



# Constants
e=1.6e-19
me=9.11e-31
# a=3e-10
hbar=1.05e-34

function createRange(start, stop, step)
    t = zeros(0)
    i = start
    while i <= stop
        append!(t, i)
        i += step
    end

    return t
end

function energyLoop(E, U, An, Bcn, step, upperLimit)
    transmissionArr = complex(zeros(0))
    reflectionArr = complex(zeros(0))
    
    @time for i in E
        grid, psi, tp, rp = nRegion(U, i, An, Bcn, step, upperLimit)
        append!(transmissionArr, tp)
        append!(reflectionArr, rp)
    end

    return transmissionArr, reflectionArr
end


function potentialStep()
    U = [0,2,0] * e
    E = 3 * e
    Bc = [1.0, 0.0]
    size_reigon = 3 
    An = [size_reigon, 5, size_reigon]

    #tarr, rarr = energyLoop(E, U, An, Bc, 0.1e-10, 0)

    grid, psi, tp, rp = nRegion(U, E[1], An, Bc, 0.01)
    plot(grid, real(psi))
    #plot!(rarr)
end



function finiteWell()
    U = [10,0,10] * e
    E = createRange(0,20000, 0.1) * e
    Bc2 = [1.0, 0.0]
    size_reigon = 3e-10 + 2e-9
    An = [size_reigon, U[1], size_reigon]

    tarr, rarr = energyLoop(E, U, An, Bc2, 1e-11, 0)
    #grid, psi, tp, rp = @time nRegion(U, 300*e, An, Bc2, 1e-11)
    
    for i in tarr
        if i == 0
            print(i)
        end
    end
    #plot!(E, rarr)

end

function equilateralTriangle()

    line_1 = createLinePotential(0.8,0)
    line_2 = createLinePotential(-0.8,0)
    Bcn = [1.0,0.0]
    
    U1, An1 = piecewisePotential(0, 2.5, 11, line_1)
    U2, An2 = piecewisePotential(0, 2.5, 11, line_2)

    display(plot(U1))
    U, An = zeros(0), zeros(0)
    U = append!(U, U1, U2)
    An = append!(An, An1, An2)

    
    
    transmissionArr = complex(zeros(0))
    reflectionArr = complex(zeros(0))
    E = 0:0.01:10
        
    @time for i in E
        grid, psi, tp, rp = nRegion(U, i, An, Bcn, 1e-3, 2)
        append!(transmissionArr, tp)
        append!(reflectionArr, rp)
    end

    plot(E, real(transmissionArr), title = "Transmission Probability", xlabel = "E", ylabel = "T(E)", label = "TP")
    #plot(E, real(reflectionArr), label = "RP")
    
end

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

function harmonicOscillator()
    new_gaussian = createHarmonicPotential()
    Bcn = [1.0, 0.0]
    U, An = piecewisePotential(-5, 5, 11, new_gaussian)

    transmissionArr = complex(zeros(0))
    reflectionArr = complex(zeros(0))
    E = 0:0.01:10
    
    #=
    @time for i in E
        grid, psi, tp, rp = nRegion(U, i, An, Bcn, 1e-2, 5)
        append!(transmissionArr, tp)
        append!(reflectionArr, rp)
    end

    plot(E, real(transmissionArr), title = "Transmission Probability", xlabel = "E", ylabel = "T(E)", label = "TP")
    #plot(E, real(reflectionArr), label = "RP")
    =#

    grid, psi, tp, rp = nRegion(U, 1, An, Bcn, 1e-2, 5)
    plotSimulation(grid, psi, 1e-7)

end

finiteWell()
