#############################################################################
#
# Functions for computing the primary radius of the bounding ellipsoid of the
# oscillatory part of the Riemann theta function.
#
# The two subroutines solve for the radius using,
#
# * Theorem 3 of [CRTF] (no derivatives)
# * A generalization of Theorems 5 and 7 of [CRTF] for N derivatives
#
# References
# ----------
#
# .. [CRTF] B. Deconinck, M.  Heil, A. Bobenko, M. van Hoeij and M. Schmies,
#    Computing Riemann Theta Functions, Mathematics of Computation, 73, (2004),
#    1417-1442.
#
# .. [DLMF] B. Deconinck, Digital Library of Mathematical Functions - Riemann
#    Theta Functions, http://dlmf.nist.gov/21
#
################################################################################

using StatsFuns, Roots

"""
    radius(epsilon::Float64,
           T::Matrix{Float64},
           derivs::Vector{Vector{Float64}} = Vector{Float64}[],
           accuracy_radius::Float64 = 5.)

Returns the primary radius of the bounding ellipsoid for computing the
Riemann theta function up to accuracy `epsilon`.

The derivative oscillatory part of the Riemann theta function has linear
growth in :math:`z` along the directions of the columns of the Riemann
matrix. `accuracy_radius` is used to determine a radius resulting in an
accurate Riemann theta for all

.. math ::

    ||z|| < \text{accuracy_radius}.

Parameters
----------
- epsilon : Requested accuracy.
- T : A gxg matrix representing the Cholesky decomposition of the imaginary
    part of a Riemann matrix.
- derivs : (Default: []) A list of directional derivatives. The number of
    directional derivatives is the order, N, of the derivative we wish to
    compute.
- accuracy_radius : (Default: 5) Radius for guaranteed region of accuracy. See above.

Returns
-------
- radius : The initial radius of the bounding ellipsoid used to truncate the
    Riemann theta function to desired accuracy.

"""
function radius(ϵ::Float64,
                T::AbstractMatrix{Float64},
                derivs::Vector{Vector{Complex128}} = Vector{Float64}[],
                accuracy_radius::Float64 = 5.)
    g = size(T,1)
    # compute the LLL-reduction of T
    U = lll_reduce(T)
    r = minimum(mapslices(norm, U, 1))

    length(derivs) == 0 && return radius0(ϵ, r, g)

    radiusN(ϵ, r, T, derivs, accuracy_radius)
end

"""
    radius0(eps::Float64, r::Float64, g::Int64)

Compute the radius with no derivatives.

Parameters
-----------
- eps : Requested accuracy.
- r : The length of the shortest lattice vector in the LLL reduction of the
    Cholesky decomposition of the imaginary part of the Riemann matrix.
- g : The genus / problem size.

Returns
-------
- radius : The initial radius of the bounding ellipsoid used to truncate the
    Riemann theta function to desired accuracy.

"""
function radius0(ϵ::Float64, r::Float64, g::Int64)
    lhs = ϵ * (2. / g) * (r / 2.)^g * gamma(g/2.)
    ins = gammainvccdf(g / 2., 1.0, lhs )
    R = sqrt(ins) + r / 2.
    S = (sqrt(2g) + r) / 2.
    max(R,S)
end


"""
   radiusN(eps, r, T, derivs, accuracy_radius=5)

Compute the radius with N derivatives.

Parameters
----------
- eps : Requested accuracy.
- r : The length of the shortest lattice vector in the LLL reduction of the
    Cholesky decomposition of the imaginary part of the Riemann matrix.
- T : A gxg matrix representing the Cholesky decomposition of the imaginary
    part of a Riemann matrix.
- derivs : A list of directional derivatives. The number of directional
    derivatives is the order, N, of the derivative we wish to compute.
- accuracy_radius : Radius for guaranteed region of accuracy. See :func:`radius`.

Returns
-------
- radius : The initial radius of the bounding ellipsoid used to truncate the
    Riemann theta function to desired accuracy.

"""
function radiusN(ϵ::Float64, r::Float64,
                 T::AbstractMatrix{Float64},
                 derivs::Vector{Vector{Complex128}},
                 accuracy_radius::Float64 = 5.)

    N, g = length(derivs), size(T,1)
    prodnormderiv = prod(norm, derivs)
    normTinv = vecnorm(inv(T))

    lhs = ϵ * r ^ g * 2 ^ (1. - g - N) /
          ( π ^ (N / 2) * g * normTinv ^ N * prodnormderiv )

    # define lower bound (guess) and attempt to solve for the radius
    lbnd = (sqrt(g + 2N + sqrt(g^2 + 8N)) + r) / 2.
    function rhs(ins::Float64)
        ai(k) = binomial(N,k) * π^(k/2.) *
                (accuracy_radius * normTinv)^k *
                gamma((g+N-k)/2.) *
                gammaccdf((g+N-k)/2., 1.0, ins)
        sum(ai, 0:N) - lhs
    end
    inszero = fzero(rhs, lbnd, 1e5)

    max(sqrt(inszero) + r / 2.0, lbnd)
end
