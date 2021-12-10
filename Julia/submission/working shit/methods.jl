

function nReigontm(E, U, boundaries)
   
    k = getTotalWaveVector(E, U)
    TM = [1 0; 0 1]
    
    for i in 1:(length(boundaries)-1)
        TM = TM * getTransferMatrix(k[i], k[i+1], boundaries[i])    
    end
    return TM    
end

function makePositive(x) 
    return x < 0 ? -x : x;
end

function modulus(z)
    return sqrt(real(z) * real(z) + imag(z) * imag(z));
end

function getTransmissionProbability(t11)
    return 1/(modulus(t11)^2)
end

function getReflectionProbability(t21, t11)
    return (modulus(t21)^2)/(modulus(t11)^2)
end

function getWaveVector(E, U)
    return sqrt(2*me*(Complex(E-U)))/hbar
end

function getTransferMatrix(k1, k2, a)

    if k1 == 0
        return [(1 - k2*a)*exp(1im *k2 *a) (1 + k2*a)*exp(-1im *k2 *a) ; (1im*k2)*exp(1im *k2 *a) (-1im*k2)*exp(-1im *k2 *a)]
    elseif k2 == 0
        return [(exp(-1im * k1 * a)/2) ((k1*a - 1im)*exp(-1im * k1 * a)/(2k1)) ; (exp(1im * k1 * a)/2) ((k1*a + 1im)*exp(1im * k1 * a)/(2k1))]    
    else
        return (1/(2*k1))*[(k1+k2)*exp(-1im*a*(k1-k2)) (k1-k2)*exp(-1im*a*(k1+k2)); (k1-k2)*exp(1im*a*(k1+k2))  (k1+k2)*exp(1im*a*(k1-k2))]
    end
end

function getWavefunction(A, B, k, x)
    return A*exp.(1im.*a.*k).+B*exp.(-1im.*a.*k)
end

function getT11(E)
    return [real(solve(i)) for i in E]
end

function getTotalWaveVector(E, U)
    return [getWaveVector(E, i) for i in U]
end