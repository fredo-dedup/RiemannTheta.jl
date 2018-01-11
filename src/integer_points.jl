###############################################################################
# Functions for computing the set of integer points used in evaluating the
# oscillatory part of the Riemann theta function.
#
# References
# ----------
#
# .. [CRTF] B. Deconinck, M.  Heil, A. Bobenko, M. van Hoeij and M. Schmies,
#    Computing Riemann Theta Functions, Mathematics of Computation, 73, (2004),
#    1417-1442.
#
# .. [DLMF] B. Deconinck, Digital Library of Mathematics Functions - Riemann Theta
#    Functions, http://dlmf.nist.gov/21
#
###############################################################################

"""
    find_int_points(g::Int64, R::Float64,
                    T::AbstractMatrix{Float64},
                    c::Vector{Float64}=zeros(g),
                    radix::Vector{Float64}=Float64[])::Vector{Vector{Float64}}

Returns the set of integer points used in computing the Riemann theta
function finite sum.

Parameters
----------
- g : Genus.
- R : Primary radius of the ellipsoid. See :func:`radius.radius` for more
information.
- T : The Cholesky decomposotion of the imaginary part of the Riemann matrix.

Returns
-------
A vector of integer vectors (given as doubles) in row-dominant
form. That is, each row of the output array is an integer vector over which
the finite sum is computed.
"""
function find_int_points(g::Int64, R::Float64,
                         T::AbstractMatrix{Float64},
                         c::Vector{Float64}=zeros(g),
                         radix::Vector{Float64}=Float64[])::Vector{Vector{Float64}}

    # determine the endpoints of this dimension of the ellipsoid and check if we
    # reached a boundary
    a = ceil(  c[g] - R / (sqrt(π) * T[g,g]) )
    b = floor( c[g] + R / (sqrt(π) * T[g,g]) )
    (a <= b) || return Vector{Float64}[]

    # construct the integer points when the final dimension is reached
    points = Vector{Float64}[]
    if g == 1
        for i in a:b  # range(int(a), int(b+1)):
            # this algorithm works backwards on the coordinates: the last
            # coordinate found is n1 if our coordinates are {n1,n2,...,ng}
            push!(points, vcat([i;] , radix))  #  5
        end
        return points

    else
        # compute new shifts, radii, radix, and recurse
        newg = g-1
        newT = T[1:end-1,1:end-1]
        newTinv = inv(newT)  #  4
        for n in a:b
            # n = first(a:b)
            chat = c[1:end-1]
            that = T[1:end-1,g]
            # newc = (chat.T - numpy.dot(newTinv, that) * (n - c[g]) ).T
            newc = chat - newTinv * that * (n - c[g])  #6
            newR = R*R - π * T[g,g]^2 * (n - c[g])^2
            newR = sqrt(newR)
            newradix = vcat([n;],radix) # 22 !!
            newpoints = _find_int_points_python(newg, newR, newT, newc, newradix)
            append!(points, newpoints)  #  3
        end
    end

    return points
end
