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

function innerpoints(T::Matrix{Float64}, radius::Float64)
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
