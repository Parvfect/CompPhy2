

using Plots

boundary = 3

function generate_angles(t, n)

    a = pi/t
    angles = []
    append!(angles, a)

    for i in 1:n
        angle = a + i*pi/(2*boundary)
        append!(angles, angle)
    end

    return angles
end

function calculate_wave_vector(angle, k)
    return sqrt((complex(-cot(angle) ^2 * k^2) + k^2))
end

function calculate_lhs(angle)
    return -cos(angle)
end

function calculate_rhs(k, K)
    return sqrt((K^2 - k^2)/k^2)
end

function new_rhs(k)
    t =[]
    for i in k:
        K = 0:0.01:10
        t_1 = [calculate_rhs(i, j) for j in K]
        append!(t, t_1)
    end
    return t
end

function new_lhs(a)

k = generate_angles(boundary, 10)
K = [calculate_wave_vector(i * boundary, i) for i in k]

K = 0:0.1:10


rhs = [calculate_rhs(i, K) for i in K]

print(k)

plot(k, real(K), label = "Solutions for Depth of a finite square well", xlabel = "k(2mE0)", ylabel = "K(2mVo)")
