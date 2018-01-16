g, num_vectors = 5, 10
z = [ rand(Complex128, g) for i in 1:num_vectors]
mode, ϵ = 0, 1e-6
accuracy_radius = 5.
tttt = rand(5,5)
Y = tttt * tttt'
Ω = Complex.(rand(g, g), Y)

############## lll_reduce.jl  ###############
b = rand(g,g)
mu = Matrix{Float64}(g,g)
B = Vector{Float64}(g)
gram_schmidt!(b, mu, B)
mu
@btime gram_schmidt!($b, $mu, $B)

############## radius.jl  ###############
@btime radius($ϵ, $(Matrix(T)), $derivs, $accuracy_radius) # 22.6 μs
@btime radius($ϵ, $T, $derivs, $accuracy_radius) # 22.6 μs


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

@btime find_int_points($g, $R, $(Matrix(T))) # 50.5 ms
@btime find_int_points($g, $R, $(T)) # 54.6 ms


T = Matrix(chol(imag.(Ω)))
n = size(T,1)
T, Rₒ, c, ns, i, nss = T, 4.5 / sqrt(π), zeros(n), Vector{Float64}(n), n, Vector{Float64}[]
@btime allrang($T, $Rₒ, $c, $ns, $i, $nss) # 6.3 μs
@time allrang(T, Rₒ, c, ns, i, nss) # 40 μs

allrang4(T, 4.5)
@btime allrang4($T, 4.5) # 4.7 μs


it = InteriorPointIterator(4.5, T)
@btime collect($it) # 118μs

n = 5
T = (tmp = rand(n,n)-0.5 ; Matrix(chol(tmp*tmp')))
Rₒ, c, ns = 5. / sqrt(π), zeros(n), Vector{Float64}(n), n
@time allrang(T, Rₒ, zeros(n), Vector{Float64}(n), n, Vector{Float64}[])
@btime allrang($T, $Rₒ, $c, $ns, $n, $(Vector{Float64}[])) # 12.4ms (37461 points)
@btime allrang4($T, 5.) # 7.2ms (37461 points)




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
