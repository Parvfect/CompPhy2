#= Solving the time Independent Schrodinger Equation for an infinite square well=#
#= No Boundary conditions on sides for momentum so starts oscillating back once the energy increases over some amount =#
using Plots
using LinearAlgebra

function get_wavefunction(n, x)
    return sqrt(2/10) * (sin((n* x * pi) / 10))
end

function get_momentum(n, k)
    return Complex( sqrt(pi * 10) * (1 - exp(- 1im * k * 10)* exp(1im * n * pi)) * n / (n^2 * pi^2 - k^2 * 10^2) )
end


function probabilityAmplitude(v)
    return v * conj(v)
end

function plot_wavefunction(n)
    grid = -5:0.02:5
    wavefunction = zeros(length(grid))
    momentum = complex(zeros(length(grid)))
    y_zero = zeros(length(grid))

    for i in 1:length(grid)
        wavefunction[i] = get_wavefunction(n, grid[i])
        momentum[i] = get_momentum(n, grid[i])
    end

    mom_abs = zeros(length(grid))
    wavefunction_abs = zeros(length(grid))

    for i in 1:length(grid)
        mom_abs[i] = probabilityAmplitude(momentum[i])
        wavefunction_abs[i] = probabilityAmplitude(wavefunction[i])
    end

    plot(grid, wavefunction_abs, xlabel = "x", label = "Wavefunction")
    plot!(grid, mom_abs, label = "Momentum")

end

plot_wavefunction(12)