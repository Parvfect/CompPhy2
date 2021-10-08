
#= Conditions , Loops, Functions
    Taking Input, Type Conversions =#

using Plots

function random_plot()
    x = 1:10
    y = rand(10)
    display(plot(x,y))
end

function hello_world_n_times(n)
    for i in 1:n
        helloWorld()
    end
end

function helloWorld()
    # Function that prints helloWorld
    print("Hello World\n")
end


random_plot()