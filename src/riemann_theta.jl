###############################################################################
#
# The primary module for computing the Riemann theta function.
#
# .. math::
#
#   Θ(z, ω) = Σ
    # \theta(z, \Omega) = \sum_{n \in \mathbb{Z}^g}
                        # e^{2 \pi i \left( \tfrac{1}{2} n \cdot \Omega n
                           # + n \cdot z \right)}
#
# References
# ----------
#
# .. [CRTF] B. Deconinck, M.  Heil, A. Bobenko, M. van Hoeij and M. Schmies,
#    Computing Riemann Theta Functions, Mathematics of Computation, 73, (2004),
#    1417-1442.
#
# .. [DLMF] B. Deconinck, Digital Library of Mathematics Functions - Riemann
#    Theta Functions, http://dlmf.nist.gov/21
#
# .. [SAGE] Computing Riemann theta functions in Sage with applications.
#    C. Swierczewski and B. Deconinck.Submitted for publication.  Available
#    online at
#    http://depts.washington.edu/bdecon/papers/pdfs/Swierczewski_Deconinck1.pdf
#
###############################################################################

"""
         oscillatory_part(z::Vector{Vector{Complex128}},
                          Ω::Matrix{Complex128},
                          mode,
                          ϵ,
                          derivs::Vector{Vector{Complex128}},
                          accuracy_radius)

Compute the oscillatory part of the Riemann theta function.


Parameters
----------
- z : A vector of complex vectors at which to evaluate the Riemann theta function.
- Omega : A Riemann matrix.
- epsilon : (Default: `1e-8`) The desired numerical accuracy.
- derivs : A vector of complex vectors giving a directional derivative.
- accuracy_radius : (Default: `5.`) The radius from the g-dimensional origin where the
    requested accuracy of the Riemann theta is guaranteed when computing
    derivatives. Not used if no derivatives of theta are requested.

Returns
-------
- The value of the Riemann theta function at each point appearing in `z`.
"""
function oscillatory_part(z::Vector{Vector{Complex128}},
                          Ω::Matrix{Complex128},
                          ϵ::Float64,
                          derivs::Vector{Vector{Complex128}},
                          accuracy_radius::Float64)
    # extract the requested information: the real part, inverse of the
    # imaginary part, and the cholesky decomposition of the imaginary part
    X = real.(Ω)
    Y = imag.(Ω)

    # In python version numpy.linalg.cholesky returns the lower triangular
    #  matrix, which is then transposed. Julia's chol returns the upper
    #  triangular matrix, hence no need to transpose.
    T = Matrix(chol(Y))

    Yinv = inv(Y)

    # compute the integer points over which we approximate the infinite sum to
    # the requested accuracy
    R = radius(ϵ, T, derivs, accuracy_radius)
    S = innerpoints(T, R)

    finite_sum(X, Yinv, T, z, S, derivs)
end


"""
         exponential_part(z::Vector{Vector{Complex128}},
                          Ω::Matrix{Complex128})

Returns the exponential part of the Riemann theta function.

This function is "vectorized" over `z`. By default, each row of `z` is
interpreted as a separate input vector to the Riemann theta function.

Parameters
----------
- z : A vector of complex vectors at which to evaluate the Riemann theta function.
- Omega : A Riemann matrix.

Returns
-------
The value of the exponential part of the Riemann theta function at
each point appearing in `z`.

"""
function exponential_part(z::Vector{Vector{Complex128}},
                          Ω::Matrix{Complex128})::Vector{Float64}
    # extract the imaginary parts of z and the inverse of the imaginary part
    # of Omega
    y = [ imag.(zz) for zz in z ]
    Yinv = inv(imag.(Ω))

    # apply the quadratic form to each vector in z
    map(y -> π * dot(y, Yinv * y), y)
end


"""
         riemanntheta(z::Vector{Vector{Complex128}},
                      Ω::Matrix{Complex128},
                      mode,
                      ϵ,
                      derivs::Vector{Vector{Complex128}},
                      accuracy_radius)::Vector{Complex128}

Returns the value of the Riemann theta function at `z` and `Omega`.

Parameters
----------
- z : A vector of complex vectors at which to evaluate the Riemann theta function.
- Omega : A Riemann matrix.
- epsilon : (Default: `1e-8`) The desired numerical accuracy.
- derivs : A vector of complex vectors giving a directional derivative.
- accuracy_radius : (Default: `5.`) The radius from the g-dimensional origin where the
    requested accuracy of the Riemann theta is guaranteed when computing
    derivatives. Not used if no derivatives of theta are requested.


Returns
-------
The value of the Riemann theta function at each point in `z`.
"""
function riemanntheta(z::Vector{Vector{Complex128}},
                      Ω::Matrix{Complex128},
                      ϵ::Float64,
                      derivs::Vector{Vector{Complex128}},
                      accuracy_radius)::Vector{Complex128}

    u = exponential_part(z, Ω)
    v = oscillatory_part(z, Ω, ϵ, derivs, accuracy_radius)

    exp.(u) .* v
end
