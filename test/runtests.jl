using RiemannTheta
using Base.Test

######### some testing values

Ω1 = [ 1.690983006 + 0.9510565162im 1.5+0.363271264im ;
      1.5+0.363271264im 1.309016994+0.9510565162im ]
T1 = Matrix(chol(imag.(Ω1)))

Ω2 = -1/(2π * im) * [ 111.207 96.616 ; 96.616 83.943 ]
T2 = Matrix(chol(imag.(Ω2)))

srand(0)
tmp = 5*rand(10,10) - 2.5
Ω3 = 5*rand(10,10) - 2.5 + (tmp * tmp') * im
T3 = Matrix(chol(imag.(Ω3)))

derivs1 = [ rand(Complex128, 2) for i in 1:1 ]
derivs2 = [ rand(Complex128, 10) for i in 1:5 ]

ϵ1, ϵ2 = 1e-3, 1e-8

##########

include("internals.jl")
include("mainfuncs.jl")
