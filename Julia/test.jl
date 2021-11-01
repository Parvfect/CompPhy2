using Test
include("wavefunction_1d.jl")
using Main.sample

@testset "Potential Step" begin

@test


#= 

Some of you who had pushed onto the later problems were finding issues with NaNs when the energy of the particle was conincidentally equal to the potential energy of the region. 

To fix this issue, you need to redefine the transfer matrix if you have this coincidence between the energy and the potential energy (E=U_i).  This could be included in your programme in the form of an "if" statement.  

The reason it arises lies in there being a different solution to the wave equation when k=0.  I have revised the session 3 worksheet (see new equations (4-6)) to explain this, so please make sure you have the latest version of the problem sheets.  Note that this won't affect your code unless E=U_i in a region of your system.
=#

end