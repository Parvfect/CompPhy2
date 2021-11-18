
include("wavefunction.jl")

function forwardPsiSim(U, E, An, Bcn, step)

    k = [getWaveVector(E, i) for i in U]
    
    len = length(k)
    # get all the boundary conditions
    Bc = complex(zeros(length(k),2))
    Bc[len,1], Bc[len,2] = Bcn[1], Bcn[2]
    boundary = An[len]
    j = len
    for i in reverse(2:length(k))
        Bc_1 = Bc[i, 1]
        Bc_2 = Bc[i, 2]
        Bcz = [Bc_1; Bc_2]
        k2 = k[i-1]
        k1 = k[i]
        t = transferMatrixMethod2(k1, k2, Bcz, boundary)
        Bc[i-1, 1], Bc[i-1,2] = t[1],t[2]
        boundary += An[i-1]
    end
    
    return k, Bc
end

function waveFunction(k, Bc, An)

    # Now depending on the resolution we have the wave vectors and we have the boundary conditions
    ptr = 0
    grid = zeros(0)
    psi = complex(zeros(0))
    for i in 1:length(k)
        
        x = ptr:0.01:(ptr+An[i])
        Bcz = [Bc[i,1], Bc[i,2]]
        y = generalisedWavefunction(x, Bcz, k[i]) 
        ptr += An[i]
        append!(grid, x)
        append!(psi, y) 
    end
    return grid, psi
end

k, Bc = forwardPsiSim([3,0,3], 0.0000001, [2,8,2], [1.0, 0.0], 1e-2)
grid2, psi = waveFunction(k, Bc, [2,8,2])
plotWavefunction(grid2, psi, 1.4)