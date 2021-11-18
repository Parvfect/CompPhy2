
function createRange(start, stop, step)
    return Array(start:step:stop)
end

function getDigits(n)
    
    t = replace(string(n), "." => "")
    
    if n < 0
        return length(t) - 1
    elseif typeof(n) == float
        return length(t) - 1
    else
        return length(t)
    end
end