################################################################################
#
# The LLL-Reduction Algorithm
#
# Authors
# -------
#
# * Grady Williams (gradyrw@gmail.com)
# * Chris Swierczewski (cswiercz@gmail.com)
#
################################################################################


"""
         gram_schmidt!(b::AbstractMatrix{Float64},
                       mu::Matrix{Float64},
                       B::Vector{Float64})

Numerically stable Gram-Schmidt algorithm.

Given a set of `n` vectors `b_1, ..., b_n` construct a collection of
orthogonal vectors `b*_1, ..., b*_n` spanning the same space.

Parameters
----------
- b : An array of `n` vectors each of size `n`.

Returns
-------
- b : Orthogonalized vectors.
- mu : Matrix of orthogonalization parameters

      mu_ij = <bi, b*j> / <b*j, b*j>
- B : Square norms of orthogonalized vectors.
"""
function gram_schmidt!(b::AbstractMatrix{Float64},
                       mu::Matrix{Float64},
                       B::Vector{Float64})
    sz = size(b,1)
    b_star = Matrix{Float64}(sz,sz)

    for i in 1:sz
        # (1) copy non-orthogonal vector: b*_i = b_i
        b_star[i,:] = b[i,:]

        # for each previously computed b*j perform the steps:
        #
        #   (2) mu_ij = <bi, b*j> / B[j]  (compute shift coefficient)
        #   (3) b*i = b*i - mu_ij*b*j     (shift by each previous b*j)
        #
        # note that this is not performed in the first iteration when i=1
        for j in 1:i-1
            # (2) compute mu_ij = <b_i, b*_j> / <b*_j, b*_j>
            numerator = dot(b[i,:], b_star[j,:])
            mu[j, i] = numerator / B[j]
            # (3) shift b*i by - mu_ij b*j
            b_star[i,:] -= mu[j,i] .* b_star[j,:]
        end

        # // (4) store the dot product Bi = <b*i, b*i>
        B[i] = dot(b_star[i,:], b_star[i,:])
    end
end


"""
  nearest_integer_shift!(mu::Matrix{Float64},
                         b::Matrix{Float64},
                         k::Int64, l::Int64, lc::Float64)

Performs operation * on p. 521 of [LLL].

Parameters
----------
- mu :
- b :
- k,l :
- lc : LLL parameter.
"""
function nearest_integer_shift!(mu::Matrix{Float64}, b::AbstractMatrix{Float64},
                                k::Int64, l::Int64, lc::Float64)
    r = round(mu[l,k])

    if (abs(mu[l,k]) > lc)
        # // shift bk by (-r*bl)
        b[k,:] -= r * b[l,:]
        # // shift mu_kj by (-r*mu_lj) for j=0,...,l-2
        mu[1:l-2,k] -= r * mu[1:l-2,l]
        # // shift mu_kl by (-r)
        mu[l,k] -= r
    end
end


"""
  lll_reduce(b₀::AbstractMatrix{Float64},
             lc::Float64=0.5, uc::Float64=0.75)

Performs Lenstra-Lenstra-Lovasv reduction on a given lattice n in
n-dimensional real space. The input matrix b is in the usual C-ordering but
the algorithm works on the columns. (This should be rewritten in a future
update for performance purposes.)

Parameters
----------
- b : Input array / `n x n` matrix.
- lc,uc : The LLL parameters.

Returns
-------
- b : The LLL reduction of the columns of the input `b`.
"""
function lll_reduce(b₀::AbstractMatrix{Float64}, lc::Float64=0.5, uc::Float64=0.75)
  b = Matrix(b₀)
  n = size(b,1)

  # initialize mu and B with zeros
  mu = zeros(n,n)
  B  = zeros(n)

  # orthogonalize the columns of b and obtain the scaling factors B and mu
  gram_schmidt!(b, mu, B)
  k = 2
  while k <= n
    nearest_integer_shift!(mu, b, k, k-1, lc)
    swap_condition = (uc - mu[k-1,k] * mu[k-1,k]) * B[k-1]

    if (B[k] < swap_condition)
      # set the "constant parameters" for this round
      mu_tmp = mu[k-1,k]
      B_tmp = B[k] + mu_tmp * mu_tmp * B[k-1]

      # scale and swap mu and B values
      mu[k-1,k] = mu_tmp * B[k-1] / B_tmp
      B[k] = B[k-1] * B[k] / B_tmp
      B[k-1] = B_tmp

      # swap b_(k-1) and b_k
      for i in 1:n
          tmp = b[k,i]
          b[k,i] = b[k-1,i]
          b[k-1,i] = tmp
      end

      # swap mu_(k-1),j and mu_k,j for j = 0,...,k-3
      for j in 1:k-3  # (j = 0; j < k-2; j++)
          tmp = mu[j,k]
          mu[j,k] = mu[j,k-1]
          mu[j,k-1] = tmp
      end

      # perform the linear transformation for i = k,...,n-1
      for i in k:n # (i = k; i < n; i++)
          tmp = mu[k-1,i] - mu_tmp * mu[k,i]
          mu[k-1,i] = mu[k,i] + mu[k-1,k] * tmp
          mu[k,i]   = tmp
      end

      if (k > 2)
          k -= 1
      end
    else
      # perform the integer shift for l = k-3, ..., 0
      for l in k-3:-1:1 #(l = k-3; l >= 0; l--)
          nearest_integer_shift!(mu, b, k, l, lc)
      end
      k += 1
    end

  end
  b
end
