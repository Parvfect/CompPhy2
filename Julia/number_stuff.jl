
using LinearAlgebra
using Plots

# println("hello world")
a = 1
dx = 0.01

x = 0:dx:1

# Brodcast Function  - needs to be looked at 
y = cos.(2pi.*x)
y1 = sin.(2pi*x)

Nx = length(x)

# println(typeof(dx))

function central_val_arr(arr)
    mid = length(arr) / 2
    if typeof(mid) == Int
        return arr[mid]
    else
        return arr[Int(mid) + 1]
    end
end

plot(x, y, label = "cosx", xlabel = "x", ylabel = 'y')

# Building a plot
p1 = plot(x, y, label = "cosx", xlabel = "x", ylabel = 'y')
p1 = plot!(y1)

# Read into the differentiate function
dydx = diff(y)./(2pi*dx)
p1 = plot!(x[1:end-1], dydx)

# Save the figure 
savefig(p1, "testplot.pdf")

display(p1)

# Handling complex numbers
z = 1im # Defines an imaginary number

print(sqrt(Complex(-1)))

#=
z = exp.[1im+2pi.*x]

p2 = plot(z, real(z), label = "Real part of z")
p2 = plot!(z, imag(z), label = "Imaginary part of z")
=#

A = [1 2 3; 4 5 6; 7 8 9]
B = rand(3,3)

C = A*B

D = A.*B

Iden = C*inv(C)
Iden - D*inv(C)

eigen(Iden)

# Defining Functions
function add(a,b)
    return a+b
end 

function add(a::String, b::String)
    return string(a, " ", b)
end

add("Hello","World")

for i in x
    println(i)
    if i>=0.5
        break
    end
end

i = 1

while i<50
    println(i)
    i+=1
end    

function suck(x::Double)
end

@time for i in x
    print(i)
end

function sum_arg(y)
    s= 0.0
    @time for i in z
        s+= i
    end
    return s
end

sum_arg(2)