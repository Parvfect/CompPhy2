using Plots
using LinearAlgebra

e, me, hbar, A3 = 1.6e-19, 9.11e-31, 1.05e-34, 1.0
save_path = "C:/Users/Parv/Documents/compphy/Julia/Data/"
U, reigon_lengths, boundaries = [0, 2, 0]*e, [2e-9, 5e-9, 2e-9], [2e-9, 7e-9]
