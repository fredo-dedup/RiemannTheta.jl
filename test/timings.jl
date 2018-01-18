using BenchmarkTools
using RiemannTheta

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

############## lll_reduce.jl  ###############

############## radius.jl  ###############

@btime radius($ϵ1, $T1) # 3.2 μs
@btime radius($ϵ1, $T1, $derivs1) # 95 μs
@btime radius($ϵ2, $T3) # 540 μs
@btime radius($ϵ2, $T3, $derivs2) # 770 μs

############## integer_points.jl  ###############

R = radius(ϵ1, T1)
@btime innerpoints(T1, R) # 4.6 μs

R = radius(ϵ2, T3)
@btime innerpoints(T3, R) # 770 μs

########### finite_sum.jl ##############

rV = 10*rand(5)
rint  = round.(rV)
rfrac = rV - rint
rint2 = round.(10*rand(5))
rX = 10*rand(5,5)
rUT = chol(rX * rX')
@btime exppart($rint, $rX, $rV, $rint2) # 295 ns


@btime normpart($rint, $rX, $rfrac) # 180ns
@btime normpart($rint, $rUT, $rfrac) # 230ns
@btime normpart2($rint, $rUT, $rfrac) # 234ns


derivs = Vector{Complex128}[]
@btime finite_sum($X, $Yinv, $(Matrix(T)), $z, $S, $derivs) # 610-690ms / 800Mb
@btime finite_sum($X, $Yinv, $T, $z, $S, $derivs) # 710-725ms / 800Mb

derivs = [ rand(Complex128, 5) for i in 1:4 ]
@btime finite_sum($X, $Yinv, $(Matrix(T)), $z, $S, $derivs) # 790ms / 800Mb
@btime finite_sum($X, $Yinv, $T, $z, $S, $derivs) # 810ms / 800Mb



##########   riemanntheta function     ###########

zs1 = [Complex128[0.5, 0.]]
@btime riemanntheta(zs1, Ω1, eps=ϵ1) # 23 μs (303 allocations: 27.34 KiB)
@btime riemanntheta(zs1, Ω1, eps=ϵ2) # 34 μs (456 allocations: 42.22 KiB)


zs2 = [ Complex128[x, 2x] for x in -1:0.01:1 ]
@btime riemanntheta(zs2, Ω1, eps=ϵ1) # 2.4 ms (31903 allocations: 2.99 MiB)
@btime riemanntheta(zs2, Ω1, eps=ϵ2) # 4.2 ms (57656 allocations: 5.39 MiB)

@btime riemanntheta(zs2, Ω1, eps=ϵ1, derivs=derivs1) # 3.4 ms (45319 allocations: 4.20 MiB)
@btime riemanntheta(zs2, Ω1, eps=ϵ2, derivs=derivs1) # 5.4 ms (71064 allocations: 6.61 MiB)


srand(0)
zs3 = [ rand(Complex128, 10) for i in 1:20 ]
@btime riemanntheta(zs3, Ω3, eps=ϵ1) # 1.1 ms (9107 allocations: 941.25 KiB)
@btime riemanntheta(zs3, Ω3, eps=ϵ2) # 3.511 ms (32507 allocations: 5.79 MiB)

@btime riemanntheta(zs3, Ω3, eps=ϵ1, derivs=derivs2) # 90.933 ms (537569 allocations: 119.46 MiB)
@btime riemanntheta(zs3, Ω3, eps=ϵ2, derivs=derivs2) # 200.545 ms (1130444 allocations: 255.16 MiB)


# comparable to openRT timings ?
# openRT at 90 ms
zs4 = [ Complex128[Complex(rand(),0.), Complex(rand(),0.)] for i in 1:10000 ]
Ω4 = [ 3im 0. ; 0. 2im ]
@btime riemanntheta(zs4, Ω4, eps=1e-8) # 167.006 ms (1420158 allocations: 133.07 MiB)
