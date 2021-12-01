#this code plots the wavefunction of an electron incident on a step potential
#it works, but it is not a good example of coding - it is uncommented - it does not use functions - it can not be easily applied to the later tasks


using Plots
using LinearAlgebra

e=1.6e-19

U1=3*e
U2=0*e
U3=3*e 
E=0.75*e

me=9.11e-31 * 0.067
a=3e-10
hbar=1.05e-34

k1=sqrt(2*me*(Complex(E-U1)))/hbar
k2=sqrt(2*me*(Complex(E-U2)))/hbar
k3=sqrt(2*me*(Complex(E-U3)))/hbar
A2=1.0

TM=(1/(2*k1))*[(k1+k2)*exp(-1im*a*(k1-k2)) (k1-k2)*exp(-1im*a*(k1+k2)); (k1-k2)*exp(1im*a*(k1+k2))  (k1+k2)*exp(1im*a*(k1-k2))]*[A2;0]
A1=TM[1]
B1=TM[2]
print(A1, " ", B1)
x=-2e-9:1e-11:a

psi1=A1*exp.(1im.*x.*k1).+B1*exp.(-1im.*x.*k1)

p1=plot(x,real(psi1))
p1=plot!(x,imag(psi1))

x=a:1e-11:2e-9

psi2=A2*exp.(1im.*x.*k2)
p1=plot!(x,real(psi2))
p1=plot!(x,imag(psi2))