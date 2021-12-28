#=
    Master Control for time independent 1D Quantum Mechancical Simulations
=#


systems = [
    "finite_well", 
    "potential_barrier",
    "gaussian_potential",
    "harmonic_oscillator",
    "double_well", 
    "quadraple_well"
]


function boundStates(system)
    if system == 2 || system == 3
        println("Your system is not compatible with calculating bound states.")
    else
        include("$(systems[system]).jl")
        plotBoundStates()
    end
end

function t11(system)
    include("$(systems[system]).jl")
    t11Loop()
end

function tprp(system)
    include("$(systems[system]).jl")
    tprpLoop()
end

function nReigon(system)

    pritnln("Enter Energy of System")
    energy = parse(Int64, readline())

    include("$(systems[system]).jl")
    nReigon(energy)
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
        println("Invalid Choice")
    end
end

main()