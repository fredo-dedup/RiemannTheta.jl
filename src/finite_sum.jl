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
#  Mods and Phase I and II functions
#  -------
#  * Stefano Carrazza, Daniel Krefl - Nov 2017
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

# Phase II special (X = 0)
function exppart_phaseII(n::Vector{Float64}, X::Matrix{Float64},
                         x::Vector{Float64})::Float64
  2π * dot(n, x)
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

function normpart2(n::Vector{Float64},
                   T::UpperTriangular{Float64,Matrix{Float64}},
                   fracshift::Vector{Float64})::Float64
    tmp1 = n + fracshift
    tmp2 = T * tmp1
    -π * dot(tmp2, tmp2)
end

#   Phase I
function normpart_phaseI(n::Vector{Float64}, T::AbstractMatrix{Float64},
                         fracshift::Vector{Float64})::Float64
    tmp1 = n + fracshift
    tmp2 = T * tmp1
    -π * dot(tmp2, tmp2)
end

#   Phase II
function normpart_phaseII(n::Vector{Float64}, T::AbstractMatrix{Float64})::Float64
    tmp2 = T * n
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


function deriv_prod_phaseI(n::Vector{Float64},
                           intshift::Vector{Float64},
                           derivs::Vector{Vector{Complex128}})::Complex128
    # compute n-intshift
    nmintshift = n - intshift

    #   Computes the dot product of each directional derivative and nmintshift.
    #   Then it computes the product of the resulting complex scalars.
    total = 1.  # real only
    for der in derivs
        term = dot(real(der), nmintshift)
        total *= term
    end

    # Compute (2*pi*i)^(nderivs) * (total_real + total_imag*i)
    pi_mult = (2π) ^ length(derivs)

    # Determines what the result of i^nderivs is, and performs the correct
    #  multiplication afterwards.
    remain = length(derivs) % 4
    (remain == 0) && return pi_mult * Complex(total, 0.)
    (remain == 1) && return pi_mult * Complex(0., total)
    (remain == 2) && return pi_mult * Complex(-total, 0.)
    (remain == 3) && return pi_mult * Complex(0., -total)
end

function deriv_prod_phaseII(n::Vector{Float64},
                            intshift::Vector{Float64},
                            derivs::Vector{Vector{Complex128}})::Complex128
    #   Computes the dot product of each directional derivative and nmintshift.
    #   Then it computes the product of the resulting complex scalars.
    total = 1.  # real only
    for der in derivs
      term = dot(real(der), n)
      total *= term
    end

    # Compute (2*pi*i)^(nderivs) * (total_real + total_imag*i)
    pi_mult = (2π) ^ length(derivs)


    # Determines what the result of i^nderivs is, and performs the correct
    #  multiplication afterwards.
    remain = length(derivs) % 4
    (remain == 0) && return pi_mult * Complex(total, 0.)
    (remain == 1) && return pi_mult * Complex(0., total)
    (remain == 2) && return pi_mult * Complex(-total, 0.)
    (remain == 3) && return pi_mult * Complex(0., -total)
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


  /*
   *   Phase I
   *   Simplified version for purely imaginary Q and z
   *
   */

  void
  finite_sum_with_derivatives_phaseI(double* fsum_real, double* fsum_imag,
				     double* X, double* Yinv, double* T,
				     double* zr, double* zi, double* S,
				     double* deriv_real, double* deriv_imag,
				     int nderivs, int g, int N, int num_vectors)
  {
    /*
      compute the shifted vectors: shift = Yinv*y as well as its integer and
      fractional parts
    */

    // ToDo: Calc first shifts then sum over points !

    // Loop over dataset
    for (int kk = 0; kk < num_vectors; kk++)
      {
	double *y = &zi[kk*g];

	int k,j;

	double intshift[g];
	double fracshift[g];
	double sum;

	for (k = 0; k < g; k++) {
          sum = 0;
          for (j = 0; j < g; j++)
	    sum += Yinv[k*g + j] * y[j];

          intshift[k] = round(sum);
          fracshift[k] = sum - intshift[k];
	}


	// compute the finite sum
	double real_total = 0, imag_total = 0;
	double npt;
	double dpr[1];
	double dpi[1];
	double* n;
	dpr[0] = 0;
	dpi[0] = 0;

	for(k = 0; k < N; k++) {
	  // the current point in S \subset ZZ^g
	  n = S + k*g;

	  // compute the "cosine" and "sine" parts of the summand

	  npt = exp(normpart(n, T, fracshift));

	  deriv_prod_phaseI(dpr, dpi, n, intshift, deriv_real, deriv_imag, nderivs, g);

	  real_total += dpr[0] * npt;
	  imag_total += dpi[0] * npt;
	}

	// store values to poiners
	fsum_real[kk] = real_total;
	fsum_imag[kk] = imag_total;
      }
  }


  /*
   *   Phase II
   *   Simplified version for purely imaginary Q and purely real z
   *
   */

  void
  finite_sum_with_derivatives_phaseII(double* fsum_real, double* fsum_imag,
				      double* X, double* Yinv, double* T,
				      double* zr, double* zi, double* S,
				      double* deriv_real, double* deriv_imag,
				      int nderivs, int g, int N, int num_vectors)
  {

    // Empty
    for (int kk = 0; kk < num_vectors; kk++) {
      fsum_real[kk] = 0;
      fsum_imag[kk] = 0;
    }

    double* n;
    double dpr[1];
    double dpi[1];

    for(int k = 0; k < N; k++) {
      // the current point in S \subset ZZ^g
      n = S + k*g;
      double npt = exp(normpart_phaseII(n, T));

      dpr[0] = 0;
      dpi[0] = 0;

      deriv_prod_phaseII(dpr, dpi, n, deriv_real, deriv_imag, nderivs, g);

      // Loop over dataset
      for (int kk = 0; kk < num_vectors; kk++)
	{
          double *x = &zr[kk*g];

          // compute the finite sum
          double ept, cpt, spt;

          // compute the "cosine" and "sine" parts of the summand
          ept = exppart_phaseII(n, X, x);
          cpt = npt*cos(ept);
          spt = npt*sin(ept);

          fsum_real[kk] += (dpr[0] * cpt - dpi[0] * spt);
          fsum_imag[kk] += (dpr[0] * spt + dpi[0] * cpt);
	}

    }
  }


  /*
   *   Phase I
   *   nth derivative over 0th derivative
   *
   */
  void
  finite_sum_with_derivatives_normalized_phaseI(double* fsum_real, double* fsum_imag,
						double* X, double* Yinv, double* T,
						double* zr, double* zi, double* S,
						double* deriv_real, double* deriv_imag,
						int nderivs, int g, int N, int num_vectors)
  {
    /*
      compute the shifted vectors: shift = Yinv*y as well as its integer and
      fractional parts
    */


    // Loop over dataset
    for (int kk = 0; kk < num_vectors; kk++)
      {
	double *y = &zi[kk*g];

	int k,j;

	double intshift[g];
	double fracshift[g];
	double sum;

	for (k = 0; k < g; k++) {
          sum = 0;
          for (j = 0; j < g; j++)
	    sum += Yinv[k*g + j] * y[j];

          intshift[k] = round(sum);
          fracshift[k] = sum - intshift[k];
	}


	// compute the finite sum
	double real_total_nom = 0, imag_total_nom = 0;
	double real_total_den = 0;

	double npt;
	double dpr[1];
	double dpi[1];
	double* n;
	dpr[0] = 0;
	dpi[0] = 0;

	for(k = 0; k < N; k++) {
	  // the current point in S \subset ZZ^g
	  n = S + k*g;

	  // compute the "cosine" and "sine" parts of the summand

	  npt = exp(normpart(n, T, fracshift));

	  deriv_prod_phaseI(dpr, dpi, n, intshift, deriv_real, deriv_imag, nderivs, g);

	  real_total_nom += dpr[0] * npt;
	  imag_total_nom += dpi[0] * npt;
	  real_total_den += npt;
	}

	fsum_real[kk] = real_total_nom/real_total_den;
	fsum_imag[kk] = imag_total_nom/real_total_den;
      }
  }


  /*
   *   Phase II
   *   Simplified version for purely imaginary Q and purely real z
   *
   */

  void
  finite_sum_with_derivatives_normalized_phaseII(double* fsum_real, double* fsum_imag,
						 double* X, double* Yinv, double* T,
						 double* zr, double* zi, double* S,
						 double* deriv_real, double* deriv_imag,
						 int nderivs, int g, int N, int num_vectors)
  {

    // Allocate temp storage
    double norm_real[num_vectors];
    double norm_imag[num_vectors];

    // Empty
    for (int kk = 0; kk < num_vectors; kk++) {
      fsum_real[kk] = 0;
      fsum_imag[kk] = 0;
      norm_real[kk] = 0;
      norm_imag[kk] = 0;
    }

    double* n;
    double dpr[1];
    double dpi[1];

    for(int k = 0; k < N; k++) {
      // the current point in S \subset ZZ^g
      n = S + k*g;
      double npt = exp(normpart_phaseII(n, T));

      dpr[0] = 0;
      dpi[0] = 0;

      deriv_prod_phaseII(dpr, dpi, n, deriv_real, deriv_imag, nderivs, g);

      // Loop over dataset
      for (int kk = 0; kk < num_vectors; kk++)
	{
          double *x = &zr[kk*g];

          // compute the finite sum
          double ept, cpt, spt;

          // compute the "cosine" and "sine" parts of the summand
          ept = exppart_phaseII(n, X, x);
          cpt = npt*cos(ept);
          spt = npt*sin(ept);

          fsum_real[kk] += (dpr[0] * cpt - dpi[0] * spt);
          fsum_imag[kk] += (dpr[0] * spt + dpi[0] * cpt);
          norm_real[kk] += cpt;
          norm_imag[kk] += spt;
	}
    }

    // Loop over dataset (setting normalization)
    for (int kk = 0; kk < num_vectors; kk++)
      {
	double old_fsum_real = fsum_real[kk];
	double norm = norm_imag[kk]*norm_imag[kk]+norm_real[kk]*norm_real[kk];

	fsum_real[kk] = (fsum_real[kk]*norm_real[kk]+fsum_imag[kk]*norm_imag[kk])/norm;
	fsum_imag[kk] = (fsum_imag[kk]*norm_real[kk]-old_fsum_real*norm_imag[kk])/norm;
      }

  }
