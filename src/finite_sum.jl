###############################################################################
#
#  Efficiently computing the finite sum part of the Riemann theta function.
#
#  Original Authors
#  -------
#  * Chris Swierczewski (@cswiercz) - September 2012, July 2016
#  * Grady Williams (@gradyrw) - October 2012
#  * Jeremy Upsal (@jupsal) - July 2016
#
###############################################################################

###############################################################################
#  exppart
#  -------
#  A helper function for the finite sum functions. Computes
#             2π < (n-intshift), (1/2)X(n-intshift) + x >
###############################################################################
function exppart(n::Vector{Float64}, X::Matrix{Float64},
                 x::Vector{Float64}, intshift::Vector{Float64})::Float64
    tmp1 = n - intshift
    tmp2 = 0.5 * X * tmp1 + x
    2π * dot(tmp1, tmp2)
end

###############################################################################
# normpart
# --------
# A helper function for the finite sum functions. Computes
#            -pi * || T*(n+fracshift) ||^2
###############################################################################
function normpart(n::Vector{Float64}, T::AbstractMatrix{Float64},
                  fracshift::Vector{Float64})::Float64
    tmp1 = n + fracshift
    tmp2 = T * tmp1
    -π * dot(tmp2, tmp2)
end



"""
      deriv_prod(n::Vector{Float64},
                 intshift::Vector{Float64},
                 derivs::Vector{Vector{Complex128}})::Complex128

Compute the real and imaginary parts of the product
         ___
         | |    2π * I <d, n-intshift>
               | |
           d in derivs
for a given n in ZZ^g.

Parameters
----------
- n : An integer vector in the finite sum ellipsoid.
- intshift : The integer part of Yinv*y.
- deriv_real, deriv_imag : The real and imaginary parts of the derivative directional vectors.

Returns
-------
- The complex "derivative product".
"""
function deriv_prod(n::Vector{Float64},
                    intshift::Vector{Float64},
                    derivs::Vector{Vector{Complex128}})::Complex128
    # compute n-intshift
    nmintshift = n - intshift

    #   Computes the dot product of each directional derivative and nmintshift.
    #   Then it computes the product of the resulting complex scalars.
    total = Complex(1., 0.)
    for der in derivs
        term = dot(der, nmintshift)
        total *= conj(term)
    end

    # Compute (2*pi*i)^(nderivs) * (total_real + total_imag*i)
    pi_mult = (2π) ^ length(derivs)

    # Determines what the result of i^nderivs is, and performs the correct
    #  multiplication afterwards.
    remain = length(derivs) % 4
    (remain == 0) && return pi_mult * total
    (remain == 1) && return pi_mult * total * Complex(0., 1.)
    (remain == 2) && return pi_mult * (-total)
    (remain == 3) && return pi_mult * total * Complex(0., -1.)
end


"""
    function finite_sum(X::Matrix{Float64},
                        Yinv::Matrix{Float64},
                        T::Matrix{Float64},
                        z::Vector{Vector{Complex128}},
                        S::Vector{Vector{Float64}},
                        derivs::Vector{Vector{Complex128}})

Computes the real and imaginary parts of the finite sum with derivatives.

Parameters
----------
X, Yinv, T : Row-major matrices such that the Riemann matrix, Omega is equal to (X +
  iY). T is the Cholesky decomposition of Y.
z : Input vectors.
S : The set of points in ZZ^g over which to compute the finite sum
derivs : the derivative directional vectors.

Returns
-------
  The finite sums for each z.
"""
function finite_sum(X::Matrix{Float64},
                    Yinv::Matrix{Float64},
                    T::AbstractMatrix{Float64},
                    z::Vector{Vector{Complex128}},
                    S::Vector{Vector{Float64}},
                    derivs::Vector{Vector{Complex128}})

    num_vectors = length(z)
    values = Vector{Complex128}(num_vectors)

    for kk in 1:num_vectors
        #   compute the shifted vectors: shift = Yinv*y as well as its integer and
        #   fractional parts
        shift = Yinv * imag(z[kk])
        intshift = round.(shift)
        fracshift = shift - intshift

        # compute the finite sum
        values[kk] = Complex(0., 0.)
        for s in S
            zpt = Complex(normpart(s, T, fracshift),
                          exppart(s, X, real(z[kk]), intshift) )
            pt = exp(zpt)
            dp = deriv_prod(s, intshift, derivs)
            values[kk] += dp * pt
        end
    end
    values
end
