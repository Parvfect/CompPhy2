#=
    Master Control for time independent 1D Quantum Mechancical Simulations
=#

include("generic_soln.jl")
include("bound_states.jl")

systems = [
    "finite_well", 
    "potential_barrier",
    "gaussian_potential",
    "harmonic_oscillator",
    "double_well", 
    "quadraple_well"
]
E = 1e-22:1e-22:3e-19


function boundStates(system)
    if system == 2 || system == 3
        println("Your system is not compatible with calculating bound states.")
    else
        include("$(systems[system]).jl")
        boundStates = getAllBoundStates(E, U, boundaries)
        print(boundStates)

    end
end

function t11(system)
    include("$(systems[system]).jl")
    t11 = energyLoopt11(E, U, boundaries)
    display(plot(real(t11)))
end

function tprp(system)
    include("$(systems[system]).jl")
    tp, rp = energyLoopTprp(E, U, boundaries)
    plot(tp)
    display(plot!(rp))
end

function nReigon(system)

    println("Enter Energy of System (will be e-21)")
    E = parse(Float64, readline()) * 10e-21
    include("$(systems[system]).jl")
    A, B = nReigon(E, U, 1.0, 0, boundaries)
    nReigonPlot(A, B, getTotalWaveVector(E, U), reigon_lengths, false)
end


function main()
    
    for (index, value) in enumerate(systems)
        println("$index. $value")
    end
    
    println("Please enter system to be simulated")
    system = parse(Int64, readline())

    println("Would you like to see the \n1. Bound States \n2.t11 variance over Energy \n3.Transmission and Reflection Probabilites \n4.Simulation at a given Energy")
    choice = parse(Int64, readline())

    if choice == 1
        boundStates(system)
    elseif choice == 2
        t11(system)
    elseif choice == 3
        tprp(system)
    elseif choice == 4
        nReigon(system)
    else
        print("Invalid Choice")
    end
    main()
end

main()