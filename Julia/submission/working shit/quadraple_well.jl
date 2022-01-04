include("generic_soln.jl")
include("bound_states.jl")

e, me, hbar, A3 = 1.6e-19, 9.11e-31, 1.05e-34, 1.0
U, reigon_lengths, boundaries = [0, -2, 0, -2, 0, -2, 0, -2, 0, -2, 0], [2e-9, 8e-9, 2e-9, 8e-9, 2e-9, 8e-9, 2e-9, 8e-9, 2e-9, 8e-9, 2e-9], [2e-9, 10e-9, 12e-9, 20e-9, 22e-9, 30e-9, 32e-9, 40e-9, 42e-9, 50e-9] 
