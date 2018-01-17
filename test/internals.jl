############## radius.jl  ###############

g = 5
U = lll_reduce(imag(Ω))
r = minimum(mapslices(norm, U, 1))

radius0(ϵ, r, g)
radiusN(ϵ, r, g, T, Vector{Complex128}[ones(g)], accuracy_radius)

radius(1e-8, T)
radius(1e-8, T, [rand(6), rand(6), rand(6)] )
radius(1e-5, T, [rand(6), rand(6), rand(6)] )
radius(1e-5, T, [rand(6)] )

radius0(eps, r, 2)
ok, égal à v Python


eps = 1e-6
U = lll_reduce(T)
r = minimum(mapslices(norm, U, 1))
g = 6
derivs = [[rand(6);]]
accuracy_radius = 5.


############## innerpoints  ###################

Ω = [ 1.690983006 + 0.9510565162im 1.5+0.363271264im ;
      1.5+0.363271264im 1.309016994+0.9510565162im ]

T = Matrix(chol(imag.(Ω)))
R = radius(T, 1e-3)
innerpoints(T, R)

ns
clipboard(ns)
[-2.0, 2.0],
[2.0, -2.0],
[0.0, 2.0],
[0.0, -2.0],
[2.0, -1.0],
[-2.0, 1.0],
[2.0, 0.0],
[-2.0, 0.0],
[0.0, -1.0],
[0.0, 1.0],
[-1.0, -1.0],
[1.0, 1.0],
[1.0, 0.0],
[-1.0, 0.0],
[-2.0, -1.0],
[2.0, 1.0],
[-1.0, 1.0],
[1.0, -1.0],
[0.0, 0.0],
[-1.0, -2.0],
[1.0, 2.0]
[-1.0, 2.0],
[1.0, -2.0],


[0.0, -2.0],
[1.0, -2.0],
[2.0, -2.0],
[-1.0, -1.0],
[0.0, -1.0],
[1.0, -1.0],
[2.0, -1.0],
[-1.0, 0.0],
[0.0, 0.0],
[1.0, 0.0],
[2.0, 0.0],
[-1.0, 1.0],
[0.0, 1.0],
[1.0, 1.0],
[2.0, 1.0],
[-1.0, 2.0],
[0.0, 2.0],
[1.0, 2.0]]


Ω = -1/(2π * im) * [ 111.207 96.616 ; 96.616 83.943 ]
T = Matrix(chol(imag.(Ω)))
R = radius(T, 1e-3)
res = innerpoints(T, R)

res = innerpoints(Matrix(chol(imag.(Ω))), 1.)

using DataVoyager
using NamedTuples
Voyager([  @NT(x=el[1], y=el[2]) for el in res ])

clipboard(res)



############## lll_reduce.jl  ###############
b = rand(5,5)
mu = Matrix{Float64}(5,5)
B = Vector{Float64}(5)
gram_schmidt!(b, mu, B)
mu
@btime gram_schmidt!($b, $mu, $B)


############## integer_points.jl  ###############
tmp = rand(3,3)
T = tmp * tmp'
T = Matrix(chol(T))
T = Matrix(T)
g = size(T, 1)
R = 2.
c = zeros(g)
radix =  Float64[]
_find_int_points_python(g, R, T, c, radix)

@btime _find_int_points_python($g, $R, $T, $c, $radix)
# g=3, R=2  (485vecs) => 370 μs
# g=5, R=5.27 (69.9k vecs) => 77 ms


########### finite_sum.jl ##############

rV = 10*rand(5)
rint  = round.(rV)
rfrac = rV - rint
rint2 = round.(10*rand(5))
rX = 10*rand(5,5)
rUT = chol(rX * rX')
@btime exppart($rint, $rX, $rV, $rint2) # 295 ns



length(derivs)
derivs = Vector{Complex128}[]
deriv_prod2(S[1000], intshift, derivs)
deriv_prod(S[1000], intshift, derivs)

derivs = [ rand(Complex128, 5) for i in 1:4 ]
deriv_prod2(S[1000], intshift, derivs)
deriv_prod(S[1000], intshift, derivs)

derivs = [ rand(Complex128, 5) for i in 1:5 ]
deriv_prod2(S[1000], intshift, derivs)
deriv_prod(S[1000], intshift, derivs)

derivs = [ rand(Complex128, 5) for i in 1:6 ]
deriv_prod2(S[1000], intshift, derivs)
deriv_prod(S[1000], intshift, derivs)

derivs = [ rand(Complex128, 5) for i in 1:7 ]
deriv_prod2(S[1000], intshift, derivs)
deriv_prod(S[1000], intshift, derivs)
