
using Pkg
Pkg.add("Plots")
using Plots


function euler_step(vec::Vector, dt, t, func::Function)
    # Takes in a vector and returns that times dt
    return dt * func(vec, t)
end

function fz(vec::Vector, t)
    
    y2 = sin(t*t) + 1/(t+1)
    return [vec[2], y2]
end



function solve(dt, t)
    y = [1.0, 0]
    y0 = []
    y1 = []

    temp = 0
    # Perform euler step
    @time while true
        y = y + euler_step(y, dt, t, fz)
        t += dt

        if temp == t
            break
        end

        append!(y0, y[1])
        append!(y1, y[2])
        display(plot(y0))
        display(plot(y0,y1))
    end
end

solve(0.01, 1000)