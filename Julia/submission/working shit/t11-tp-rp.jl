include("methods.jl")

function nReigont11(E, U, boundaries)
   
    k = getTotalWaveVector(E, U)
    TM = [1 0; 0 1]
    
    for i in 1:(length(boundaries)-1)
        TM = TM * getTransferMatrix(k[i], k[i+1], boundaries[i])    
    end
    return TM    
end