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
- mode: 0 no restrictions on z, Ω
        1 phaseI  : Ω imaginary, z imaginary
        2 phaseII : Ω imaginary, z real
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
                          mode,
                          ϵ,
                          derivs::Vector{Vector{Complex128}},
                          accuracy_radius)
    # extract the requested information: the real part, inverse of the
    # imaginary part, and the cholesky decomposition of the imaginary part
    X = real(Ω)
    Y = imag(Ω)

    # In python version numpy.linalg.cholesky returns the lower triangular
    #  matrix, which is then transposed. Julia's chol returns the upper
    #  triangular matrix, hence no need to transpose.
    T = Matrix(chol(Y))

    Yinv = mode != 2 ? inv(Y) : zeros(Y)

    # compute the integer points over which we approximate the infinite sum to
    # the requested accuracy
    R = radius(ϵ, T, derivs, accuracy_radius)
    S = find_int_points(g, R, T)

    # compute the finite sum for each z-vector
    if mode == 0
        values = finite_sum(X, Yinv, T, z, S, derivs)
    elseif mode == 1
        values = finite_sum_phaseI(X, Yinv, T, z, S, derivs)
    elseif mode == 2
        values = finite_sum_phaseII(X, Yinv, T, z, S, derivs)
    end

    values
end


"""
         exponential_part(z::Vector{Vector{Complex128}},
                          Ω::Matrix{Complex128})

Returns the exponential part of the Riemann theta function.

This function is "vectorized" over `z`. By default, each row of `z` is
interpreted as a separate input vector to the Riemann theta function.

Parameters
----------
- z : A complex row-vector or list of row-vectors at which to evaluate the
Riemann theta function.
- Omega : A Riemann matrix.

Returns
-------
The value of the exponential part of the Riemann theta function at
each point appearing in `z`.

"""
function exponential_part(z::Vector{Vector{Complex128}},
                          Ω::Matrix{Complex128})
    # extract the imaginary parts of z and the inverse of the imaginary part
    # of Omega
    y = imag.(z)
    Yinv = inv(imag.(Ω))

    # apply the quadratic form to each vector in z
    # quad = lambda yi: numpy.pi*numpy.dot(yi,numpy.dot(Yinv,yi))
    # exponents = numpy.apply_along_axis(quad, axis, y)
end


"""
         riemanntheta(z::Vector{Vector{Complex128}},
                      Ω::Matrix{Complex128},
                      mode,
                      ϵ,
                      derivs::Vector{Vector{Complex128}},
                      accuracy_radius)

Returns the value of the Riemann theta function at `z` and `Omega`.

In many applications it's preferred to use :meth:`exponential_part` and
:meth:`oscillatory_part` due to the double-exponential growth of theta
in the directions of the columns of `Omega`.

Parameters
----------
z : complex[:]
A complex row-vector or list of row-vectors at which to evaulate the
Riemann theta function.
Omega : complex[:,:]
A Riemann matrix.

Returns
-------
The value of the Riemann theta function at each point in `z`.
"""
function riemanntheta(z::Vector{Vector{Complex128}},
                      Ω::Matrix{Complex128},
                      mode,
                      ϵ,
                      derivs::Vector{Vector{Complex128}},
                      accuracy_radius)

    u = exponential_part(z, Ω)
    v = oscillatory_part(z, Ω, mode, ϵ, derivs, accuracy_radius)

    exp.(u) .* v
end



######### normalized_eval

"""

       normalized_oscillatory_part(z::Vector{Vector{Complex128}},
                                   Ω::Matrix{Complex128},
                                   mode,
                                   ϵ,
                                   derivs::Vector{Vector{Complex128}},
                                   accuracy_radius)

Compute the normalized oscillatory part of the Riemann theta function.

The oscillatory part of the Riemann theta function is the infinite
summation left over after factoring out the double-exponential growth.

This function is "vectorized" over `z` in order to take advantage of the
uniform approximation theorem.

Parameters
----------
- z : A vector of complex vectors at which to evaluate the Riemann theta function.
- Omega : A Riemann matrix.
- mode: 0 no restrictions on z, Ω
      1 phaseI  : Ω imaginary, z imaginary
      2 phaseII : Ω imaginary, z real
- epsilon : (Default: `1e-8`) The desired numerical accuracy.
- derivs : A vector of complex vectors giving a directional derivative.
- accuracy_radius : (Default: `5.`) The radius from the g-dimensional origin where the
  requested accuracy of the Riemann theta is guaranteed when computing
  derivatives. Not used if no derivatives of theta are requested.

Returns
-------
- The value of the Riemann theta function at each point appearing in `z`.

"""
# function normalized_oscillatory_part(z::Vector{Vector{Complex128}},
#                             Ω::Matrix{Complex128},
#                             mode,
#                             ϵ,
#                             derivs::Vector{Vector{Complex128}},
#                             accuracy_radius)
#
#
#     # coerce z to a numpy array and determine the problem size: the genus g and
#     # the number of vectors to compute
#     z = numpy.array(z, dtype=numpy.complex)
#     if len(z.shape) == 1:
#         num_vectors = 1
#         g = len(z)
#     else:
#         num_vectors = z.shape[abs(axis-1)]
#         g = z.shape[axis]
#     z = z.flatten()
#
#     # coerce Omega to a numpy array and extract the requested information: the
#     # real part, inverse of the imaginary part, and the cholesky decomposition
#     # of the imaginary part
#     Omega = numpy.array(Omega, dtype=numpy.complex)
#     Y = Omega.imag
#     X = numpy.ascontiguousarray(Omega.real)
#
#     if(Y.shape[0]!=1):
#         _T = numpy.linalg.cholesky(Y).T
#     else:
#         _T = numpy.sqrt(Y)
#
#     T = numpy.ascontiguousarray(_T)
#
#     if(mode!=2):
#         Yinv = numpy.ascontiguousarray(numpy.linalg.inv(Y))
#     else:
#         Yinv = numpy.ascontiguousarray(numpy.zeros(Y.shape))
#
#
#     # compute the integer points over which we approximate the infinite sum to
#     # the requested accuracy
#     if isinstance(derivs, list):
#
#         R = radius(epsilon, _T, derivs=derivs, accuracy_radius=accuracy_radius)
#         S = numpy.ascontiguousarray(integer_points_python(g,R,_T))
#         N = S.shape[0]
#
#         # set up storage locations and vectors
#         real = <double*>malloc(sizeof(double)*num_vectors)
#         imag = <double*>malloc(sizeof(double)*num_vectors)
#
#         values = numpy.zeros(num_vectors, dtype=numpy.complex)
#
#         x = numpy.ascontiguousarray(z.real, dtype=numpy.double)
#         y = numpy.ascontiguousarray(z.imag, dtype=numpy.double)
#
#         # get the derivatives
#         if len(derivs):
#             derivs = numpy.array(derivs, dtype=numpy.complex).flatten()
#             nderivs = len(derivs) / g
#             derivs_real = numpy.ascontiguousarray(derivs.real, dtype=numpy.double)
#             derivs_imag = numpy.ascontiguousarray(derivs.imag, dtype=numpy.double)
#
#             # compute the finite sum for each z-vector
#             if(mode==0):
#                 values_nom = numpy.zeros(num_vectors, dtype=numpy.complex)
#                 values_den = numpy.zeros(num_vectors, dtype=numpy.complex)
#
#                 finite_sum_with_derivatives(real, imag,
#                                             &X[0,0], &Yinv[0,0], &T[0,0],&x[0],
#                                             &y[0], &S[0,0],
#                                             &derivs_real[0], &derivs_imag[0],nderivs, g, N, num_vectors)
#
#                 for k in range(num_vectors):
#                     values_nom[k] = numpy.complex(real[k] + 1.0j*imag[k])
#
#                 finite_sum_without_derivatives(real, imag,
#                                                &x[0], &y[0], &X[0,0],&Yinv[0,0], &T[0,0], &S[0,0],g, N, num_vectors)
#                 for k in range(num_vectors):
#                     values_den[k] = numpy.complex(real[k] + 1.0j*imag[k])
#
#                 values = values_nom/values_den
#
#             elif(mode==1):
#                 finite_sum_with_derivatives_normalized_phaseI(real, imag,
#                                             &X[0,0], &Yinv[0,0], &T[0,0],&x[0],
#                                             &y[0], &S[0,0],
#                                             &derivs_real[0], &derivs_imag[0],nderivs, g, N, num_vectors)
#
#                 for k in range(num_vectors):
#                     values[k] = numpy.complex(real[k] + 1.0j*imag[k])
#
#             elif(mode==2):
#                 finite_sum_with_derivatives_normalized_phaseII(real, imag,
#                                             &X[0,0], &Yinv[0,0], &T[0,0],&x[0],
#                                             &y[0], &S[0,0],
#                                             &derivs_real[0], &derivs_imag[0],nderivs, g, N, num_vectors)
#
#                 for k in range(num_vectors):
#                     values[k] = numpy.complex(real[k] + 1.0j*imag[k])
#
#         else:
#             return numpy.ones(z.shape);
#
#
#         free(real)
#         free(imag)
#
#         if num_vectors == 1:
#             return values[0]
#         else:
#             return values
#
#     else:
#
#         output = numpy.zeros(shape=(derivs.shape[0], num_vectors), dtype=numpy.complex)
#
#         # compute the finite sum for each z-vector
#         if mode == 0:
#
#             real = <double*>malloc(sizeof(double)*num_vectors)
#             imag = <double*>malloc(sizeof(double)*num_vectors)
#
#             for i, value in enumerate(derivs):
#                 # get the derivatives
#                 if len(value):
#                     value = numpy.array(value, dtype=numpy.complex).flatten()
#                     nderivs = len(value) / g
#                     derivs_real = numpy.ascontiguousarray(value.real, dtype=numpy.double)
#                     derivs_imag = numpy.ascontiguousarray(value.imag, dtype=numpy.double)
#
#                     R = radius(epsilon, _T, derivs=value, accuracy_radius=accuracy_radius)
#                     S = numpy.ascontiguousarray(integer_points_python(g,R,_T))
#                     N = S.shape[0]
#
#                     # set up storage locations and vectors
#                     values = numpy.zeros(num_vectors, dtype=numpy.complex)
#
#                     x = numpy.ascontiguousarray(z.real, dtype=numpy.double)
#                     y = numpy.ascontiguousarray(z.imag, dtype=numpy.double)
#
#                     values_nom = numpy.zeros(num_vectors, dtype=numpy.complex)
#                     values_den = numpy.zeros(num_vectors, dtype=numpy.complex)
#
#                     finite_sum_with_derivatives(real, imag,
#                                                 &X[0,0], &Yinv[0,0], &T[0,0],&x[0],
#                                                 &y[0], &S[0,0],
#                                                 &derivs_real[0], &derivs_imag[0],nderivs, g, N, num_vectors)
#
#                     for k in range(num_vectors):
#                         values_nom[k] = numpy.complex(real[k] + 1.0j*imag[k])
#
#                     finite_sum_without_derivatives(real, imag,
#                                                    &x[0], &y[0], &X[0,0],&Yinv[0,0], &T[0,0], &S[0,0],g, N, num_vectors)
#                     for k in range(num_vectors):
#                         values_den[k] = numpy.complex(real[k] + 1.0j*imag[k])
#
#                     values = values_nom/values_den
#
#                 else:
#                     free(real)
#                     free(imag)
#                     return numpy.ones(z.shape)
#
#                 if num_vectors == 1:
#                     output[i] = values[0]
#                 else:
#                     output[i] = values
#
#             free(real)
#             free(imag)
#         else:
#
#             R = radius(epsilon, _T, derivs=max(derivs, key=len), accuracy_radius=accuracy_radius)
#             S = numpy.ascontiguousarray(integer_points_python(g,R,_T))
#             N = S.shape[0]
#
#             # set up storage locations and vectors
#             n_derivs = numpy.ascontiguousarray([len(item) for item in derivs], dtype=numpy.int32)
#             inderivs = numpy.array([item for sublist in derivs for item in sublist], dtype=numpy.complex)
#             derivs_real = numpy.ascontiguousarray(inderivs.real, dtype=numpy.double)
#             derivs_imag = numpy.ascontiguousarray(inderivs.imag, dtype=numpy.double)
#
#             real = <double*>malloc(sizeof(double)*num_vectors*derivs.shape[0])
#             imag = <double*>malloc(sizeof(double)*num_vectors*derivs.shape[0])
#
#             if mode == 1:
#                 y = numpy.ascontiguousarray(z.imag, dtype=numpy.double)
#                 finite_sum_with_multi_derivatives_normalized_phaseI(real, imag, &Yinv[0,0], &T[0,0], &y[0], &S[0,0],
#                                                                     &derivs_real[0], &derivs_imag[0], &n_derivs[0],
#                                                                     len(derivs), g, N, num_vectors)
#             if mode == 2:
#                 x = numpy.ascontiguousarray(z.real, dtype=numpy.double)
#                 finite_sum_with_multi_derivatives_normalized_phaseII(real, imag, &X[0,0], &T[0,0], &x[0], &S[0,0],
#                                                                      &derivs_real[0], &derivs_imag[0], &n_derivs[0],
#                                                                      len(derivs), g, N, num_vectors)
#
#             for i in range(derivs.shape[0]):
#                     if num_vectors == 1:
#                         output[i] = numpy.complex(real[i] + 1.0j*imag[i])
#                     else:
#                         for k in range(num_vectors):
#                             output[i, k] = numpy.complex(real[k+i*num_vectors] + 1.0j*imag[k+i*num_vectors])
#
#             free(real)
#             free(imag)
#
#         return output
# end

    #
    # def oscillatory_part_gradient(self, z, Omega, mode=0, epsilon=1e-8,
    #                               accuracy_radius=5, axis=1, **kwds):
    #     r"""Returns the oscillatory part of the gradient of Riemann theta.
    #
    #     A helper function for :meth:`gradient`. Useful in it's own right for
    #     detecting if the gradient vanishes since the exponential part of
    #     Riemann theta never does as well as for controlling exponential growth in
    #     rational functions of theta.
    #
    #     Parameters
    #     ----------
    #     (see :meth:`oscillatory_part`)
    #
    #     Returns
    #     -------
    #     array
    #         If a single z-argument is given then returns a 1-dimensional Numpy
    #         array representing the gradient. If multiple z-arguments are given
    #         then returns a 2-dimensional Numpy array of gradients where each
    #         gradient is listed along `axis`.
    #     """
    #     # get the genus and number of z-vectors
    #     Omega = numpy.array(Omega, dtype=numpy.complex)
    #     g = Omega.shape[0]
    #     gradients = numpy.zeros_like(z, dtype=numpy.complex)
    #
    #     for i in range(g):
    #         # construct the direction derivative \partial z_i and compute the
    #         # derivatives in that direction
    #         derivs = [[0]*g]
    #         derivs[0][i] = 1
    #         partial_zi = oscillatory_part(z, Omega, mode, epsilon, derivs,
    #                                       accuracy_radius, axis)
    #
    #         # if axis=1 then store gradients in rows (since input z-vectors are
    #         # given as rows of a matrix). otherwise, store column-wise
    #         if axis:
    #             gradients[:,i] = partial_zi
    #         else:
    #             gradients[i,:] = partial_zi
    #
    #     return gradients
    #
    # def gradient(self, z, Omega, mode=0, epsilon=1e-8, accuracy_radius=5,
    #              axis=1, **kwds):
    #     r"""Returns the gradient of Riemann theta.
    #
    #     The gradient of :math:`\theta((z_1, \ldots, z_g), \Omega)` is defined
    #     to be
    #
    #     .. math::
    #
    #         \grad \theta(z,\Omega) = \left( \ldots, \partial_{z_i} \theta(z,\Omega), \ldots \right)
    #
    #     where :math:`\partial_{z_i}` denotes the directional derivative in the
    #     :math:`[\ldots,0,1,0,\ldots]` direction. (The "1" occurs in the
    #     :math:`i` th position.)
    #
    #     See Also
    #     --------
    #     oscillatory_part_gradient
    #
    #     Parameters
    #     ----------
    #     (see :meth:`oscillatory_part`)
    #
    #     Returns
    #     -------
    #     array
    #         If a single z-argument is given then returns a 1-dimensional Numpy
    #         array representing the gradient. If multiple z-arguments are given
    #         then returns a 2-dimensional Numpy array of gradients where each
    #         gradient is listed along `axis`.
    #
    #     Notes
    #     -----
    #     This could be made shorter by using :meth:`eval` but we choose not to
    #     for performance reasons. (We avoid having to compute the exponential
    #     part multiple times.)
    #
    #     """
    #     # compute the gradient of the oscillatory part, the exponential part,
    #     # and combine "axis-wise": if the z-vectors are given as row vectors
    #     # then multiply each column of the gradient of the oscillatory part by
    #     # the exponential part
    #     osc_gradients = self.oscillatory_part_gradient(
    #         z, Omega, mode,epsilon, accuracy_radius, axis, **kwds)
    #     exp_parts = self.exponential_part(z, Omega, **kwds)
    #     exp_parts = numpy.exp(exp_parts)
    #     gradients = numpy.apply_along_axis(
    #         lambda dz: exp_parts*dz,
    #         (axis+1)%2,
    #         osc_gradients)
    #     return gradients
    #
    # def oscillatory_part_hessian(self, z, Omega,mode=0, epsilon=1e-8,
    #                              accuracy_radius=5, axis=1, **kwds):
    #     r"""Returns the oscillatory part of the Hessian of Riemann theta.
    #
    #     A helper function for :meth:`hessian`.  Useful in it's own right for
    #     detecting if the Hessian vanishes since the exponential part of Riemann
    #     theta never does as well as for controlling exponential growth in
    #     rational functions of theta.
    #
    #     Parameters
    #     ----------
    #     (see :meth:`oscillatory_part`)
    #
    #     Returns
    #     -------
    #     array
    #         If a single z-argument is givne then returns a 2-dimensional Numpy
    #         array representing the Hessian. If multiple z-arguments are given
    #         then returns a 3-dimensional Numpy array of Hessians where each
    #         Hessian is indexed by the 0th coordinate.
    #     """
    #     # get the genus and number of z-vectors
    #     Omega = numpy.array(Omega, dtype=numpy.complex)
    #     g = Omega.shape[0]
    #     z = numpy.array(z, dtype=numpy.complex)
    #     if len(z.shape) == 1:
    #         n = 1
    #     else:
    #         n = z.shape[(axis+1) % 2]
    #     hessians = numpy.zeros((n,g,g), dtype=numpy.complex)
    #
    #     # since Riemann theta is analytic we only need to compute the lower
    #     # triangular portion of the Hessian and symmetrize
    #     lower_derivs = []
    #     for i in range(g):
    #         d1 = [0]*g
    #         d1[i] = 1
    #         for j in range(i):
    #             d2 = [0]*g
    #             d2[j] = 1
    #
    #             # construct the derivative and evaluate the corresponding zi zj
    #             # derivative across all input vectors
    #             derivs = [d1, d2]
    #             partial_zizj = oscillatory_part(z, Omega, mode, epsilon, derivs,
    #                                             accuracy_radius, axis)
    #             hessians[:,i,j] = partial_zizj
    #
    #     # symmetrize
    #     hessians += hessians.transpose((0,2,1))
    #
    #     # finally, compute the diagonal entries of each hessian
    #     for i in range(g):
    #         d1 = [0]*g
    #         d1[i] = 1
    #         derivs = [d1,d1]
    #         partial_zizi = oscillatory_part(z, Omega, mode, epsilon, derivs,
    #                                         accuracy_radius, axis, mode)
    #         hessians[:,i,i] = partial_zizi
    #
    #     if n == 1:
    #         return hessians[0,:,:]
    #     else:
    #         return hessians
