include("generic_soln.jl")
include("bound_states.jl")

e, me, hbar, A3 = 1.6e-19, 9.11e-31, 1.05e-34, 1.0
save_path = "C:/Users/Parv/Documents/compphy/Julia/Data/"
U, reigon_lengths, boundaries = [0, -3, 0]*e, [2e-9, 8e-9, 2e-9], [2e-9, 10e-9]



function nReigonSim(E)
    A, B = nReigon(E, U, 1.0, 0, boundaries)
    nReigonPlot(A, B, getTotalWaveVector(E, U), reigon_lengths, false)
end

function t11Sim()
    E = 1e-22:1e-22:3e-19
    energyLoop(E, U, boundaries)
end

function plotBoundStates()
    E = 1e-22:1e-22:3e-19
    boundStates = getAllBoundStates(E, U, boundaries)
    print(boundStates)
    nReigonSim(boundStates[6])
end

# Bound States - Highest to Lowest
#E = 3.843359e-19
#E = 1.2678747928e-19
#E = 7.1410000239e-20
#E = 1.4122357142399997e-20
#E = 3.531310235e-21
#E = 0.8828720751e-21

#nReigonSim(3.5e-21)
#t11Sim()
E = 1e-28:1e-22:5e-19
#t11 = energyLoop(E, U, boundaries)
#display(plot(real(t11)))
plotBoundStates()