################################################################################
#
#   Iterator for the integer coordinates in ℤⁿ lying inside the n
#          dimensional ellipsoid
#
#   Implementation of Remarks page 1426 of
#   B. Deconinck, M.  Heil, A. Bobenko, M. van Hoeij and M. Schmies,
#    Computing Riemann Theta Functions, Mathematics of Computation, 73, (2004),
#    1417-1442.
#
################################################################################

using BenchmarkTools

struct InteriorPointIterator
    radius::Float64
    T::Matrix{Float64}
end


################################################################################

T = (tmp = 10*rand(n,n)-5. ; Matrix(chol(tmp*tmp')))

Rₒ, c, ns, i = 4.5 / sqrt(π), zeros(n), Vector{Float64}(n), n
function allrang(T, Rₒ, c, ns, i, nss)
    hw = Rₒ / T[i,i]
    δc = inv(T[1:i-1,1:i-1]) * T[1:i-1,i]
    for ng in ceil(c[i]-hw):floor(c[i]+hw)
        # ng = first(ceil(c[i]-hw):floor(c[i]+hw))
        ns[i] = ng
        if i == 1
            push!(nss, copy(ns))
        else
            δcn = (ng - c[i])
            nc = c[1:i-1] - δc * δcn
            nRₒ = sqrt( Rₒ^2 - (T[i,i] * δcn)^2 )
            allrang(T, nRₒ, nc, ns, i-1, nss)
        end
    end
    nss
end

ns = Vector{Float64}[]
allrang(T, 4.5 / sqrt(π), zeros(n), Vector{Float64}(n), n, Vector{Float64}[])
ns

Ω = [ 1.690983006 + 0.9510565162im 1.5+0.363271264im ;
      1.5+0.363271264im 1.309016994+0.9510565162im ]


T = Matrix(chol(imag.(Ω)))


n = size(Ω,1)
@btime allrang(T, 4.5 / sqrt(π), zeros(n), Vector{Float64}(n), n, Vector{Float64}[]) # 7-8 μs
ns
clipboard(ns)

typeof(ceil(c[i]-hw):floor(c[i]+hw))

####################

Ω = [ 1.690983006 + 0.9510565162im 1.5+0.363271264im ;
      1.5+0.363271264im 1.309016994+0.9510565162im ]

ipi = InteriorPointIterator(4.5, Matrix(chol(imag.(Ω))))

function update_from_i!(rits, tggs, δcs, Rₒs, cs, ns, i)
    hw   = Rₒs[i] / tggs[i]
    iter = ceil(cs[i][i]-hw):floor(cs[i][i]+hw)
    state = start(iter)
    if i==1
        rits[i] = (iter, state)
    else # recurse if not fully updated
        ns[i], state = next(iter, state)
        rits[i] = (iter, state)
        δcn = (ns[i] - cs[i][i])
        cs[i-1] = cs[i][1:i-1] - δcs[i] * δcn
        Rₒs[i-1] = sqrt( Rₒs[i]^2 - (tggs[i] * δcn)^2 )
        update_from_i!(rits, tggs, δcs, Rₒs, cs, ns, i-1)
    end
end

function next_to_i(rits, tggs, δcs, Rₒs, cs, ns, i)
    if done(rits[i]...)
        next_to_i(rits, tggs, δcs, Rₒs, cs, ns, i+1)
    else
        ns[i], state = next(rits[i]...)
        rits[i] = (rits[i][1], state)
        δcn = (ns[i] - cs[i][i])
        cs[i-1] = cs[i][1:i-1] - δcs[i] * δcn
        Rₒs[i-1] = sqrt( Rₒs[i]^2 - (tggs[i] * δcn)^2 )
        update_from_i!(rits, tggs, δcs, Rₒs, cs, ns, i-1)
    end
end


function Base.start(ipi::InteriorPointIterator)
    n    = size(ipi.T, 1)
    δcs  = [ inv(ipi.T[1:i-1,1:i-1]) * ipi.T[1:i-1,i] for i in 1:n ]
    tggs = diag(ipi.T)
    rits = Vector{Tuple{Range{Float64}, Int64}}(n) # range iterators
    Rₒs  = Vector{Float64}(n)
    Rₒs[n] = ipi.radius / sqrt(π)
    cs    = Vector{Vector{Float64}}(n)
    cs[n] = zeros(n)
    ns    = Vector{Float64}(n)

    update_from_i!(rits, tggs, δcs, Rₒs, cs, ns, n)

    (rits, tggs, δcs, Rₒs, cs, ns)
end

function Base.next(ipi::InteriorPointIterator, state)
    rits, tggs, δcs, Rₒs, cs, ns = state

    done(rits[1]...) && next_to_i(rits, tggs, δcs, Rₒs, cs, ns, 2)

    ns[1], istate = next(rits[1]...)
    rits[1] = (rits[1][1], istate)
    (copy(ns), state)
end

function Base.done(ipi::InteriorPointIterator, state)
    all(t -> done(t...), state[1])
end

Base.eltype(::Type{InteriorPointIterator}) = Vector{Float64}
Base.iteratorsize(::Type{InteriorPointIterator}) = Base.SizeUnknown()

it = InteriorPointIterator(4.5, Matrix(chol(imag.(Ω))))
st = start(it)
ns, st = next(it, st)
ns
ns, st = next(it, st)
ns
ns, st = next(it, st)
ns
ns, st = next(it, st)
ns
ns, st = next(it, st)
ns



@btime res = collect(InteriorPointIterator(4.5, Matrix(chol(imag.(Ω))))) # 130 μs

@btime (for ns in InteriorPointIterator(4.5, Matrix(chol(imag.(Ω)))) ; end)



################# timings  ########################

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


it = InteriorPointIterator(4.5, T)
@btime collect($it) # 118μs






@btime collect(InteriorPointIterator(4.5, T))

ns = Vector{Float64}[]
@btime allrang(T, 4.5 / sqrt(π), zeros(n), Vector{Float64}(n), n)



using DataVoyager
using NamedTuples
Voyager([  @NT(x=el[1], y=el[2]) for el in res ])

clipboard(res)


################################################################################

function allrang2(T, Rₒ, c, ns, i, chan::Channel)
    hw = Rₒ / T[i,i]
    δc = inv(T[1:i-1,1:i-1]) * T[1:i-1,i]
    for ng in ceil(c[i]-hw):floor(c[i]+hw)
        # ng = first(ceil(c[i]-hw):floor(c[i]+hw))
        ns[i] = ng
        if i == 1
            put!(chan, copy(ns))
        else
            δcn = (ng - c[i])
            nc = c[1:i-1] - δc * δcn
            nRₒ = sqrt( Rₒ^2 - (T[i,i] * δcn)^2 )
            allrang2(T, nRₒ, nc, ns, i-1, chan)
        end
    end
end

function allrang3(T, Rₒ, c, ns, i)
    chan::Channel -> allrang2(T, Rₒ, c, ns, i, chan)
end


chn = Channel(allrang3(T, Rₒ, c, ns, i))
take!(chn)
take!(chn)
take!(chn)
take!(chn)
take!(chn)
take!(chn)
take!(chn)
take!(chn)
collect(chn)

@btime collect(chn)


T = Matrix(chol(imag.(Ω)))
n = size(T,1)
@btime allrang(T, 4.5 / sqrt(π), zeros(n), Vector{Float64}(n), n, Vector{Float64}[]) # 7 μs
T, Rₒ, c, ns, i, nss = T, 4.5 / sqrt(π), zeros(n), Vector{Float64}(n), n, Vector{Float64}[]
chn = Channel(allrang3(T, Rₒ, c, ns, i))
collect(chn)
@btime collect($chn) # 6.3 μs
@btime collect($(Channel(allrang3(T, Rₒ, c, ns, i)))) # 6.3 μs

@time collect(Channel(allrang3(T, Rₒ, c, ns, i)))


############################################################################



function _innerpoints(Tgg, δcs, Rₒ, c, ns, g)
    hw = Rₒ / Tgg[g]
    for ng in ceil(c[g]-hw):floor(c[g]+hw)
        ns[g] = ng
        if i == 1
            push!(nss, copy(ns))
        else
            δcn = (ng - c[g])
            nc = c[1:i-1] - δcs[g] * δcn
            nRₒ = sqrt( Rₒ^2 - (Tgg[g] * δcn)^2 )
            allrang(T, nRₒ, nc, ns, i-1, nss)
        end
    end
    nss
end

function allrang4(T::Matrix{Float64}, radius::Float64)
    n   = size(T, 1)
    δcs = Vector{Float64}[ inv(T[1:i-1,1:i-1]) * T[1:i-1,i] for i in 1:n ]
    Tgg = diag(T)
    ns     = Vector{Float64}(n)
    points = Vector{Float64}[]

    function _innerpoints(Rₒ, c, g)
        hw = Rₒ / Tgg[g]
        for ng in ceil(c[g]-hw):floor(c[g]+hw)
            ns[g] = ng
            if g == 1
                push!(points, copy(ns))
            else
                δcn = (ng - c[g])
                nc = c[1:g-1] - δcs[g] * δcn
                nRₒ = sqrt( Rₒ^2 - (Tgg[g] * δcn)^2 )
                _innerpoints(nRₒ, nc, g-1)
            end
        end
    end

    _innerpoints(radius / sqrt(π), zeros(n), n)

    points
end


allrang4()
