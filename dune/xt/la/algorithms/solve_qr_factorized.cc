// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)
//   Tobias Leibner  (2017)

#include "config.h"

#include <algorithm>
#include <cmath>
#include <iostream>

#include <dune/common/exceptions.hh>

#include "solve_sym_tridiag_posdef.hh"

#include <dune/xt/common/lapacke.hh>

namespace Dune {
namespace XT {
namespace LA {

// Solves Ax = b, where A is lower triangular. b is overwritten with the solution x.
void solve_lower_triangular(const double* A, double* b, int rows)
{
  Common::Blas::dtrsm(Common::Blas::row_major(),
                      Common::Blas::left(),
                      Common::Blas::lower(),
                      Common::Blas::no_trans(),
                      Common::Blas::non_unit(),
                      rows,
                      1,
                      1.,
                      A,
                      rows,
                      b,
                      1);
}

// Solves A^T x = b, where A is lower triangular. b is overwritten with the solution x.
void solve_lower_triangular_transposed(const double* A, double* b, int rows)
{
  Common::Blas::dtrsm(Common::Blas::row_major(),
                      Common::Blas::left(),
                      Common::Blas::lower(),
                      Common::Blas::trans(),
                      Common::Blas::non_unit(),
                      rows,
                      1,
                      1.,
                      A,
                      rows,
                      b,
                      1);
}

// Solves Ax = b, where A is upper triangular. b is overwritten with the solution x.
void solve_upper_triangular(const double* A, double* b, int rows)
{
  Common::Blas::dtrsm(Common::Blas::row_major(),
                      Common::Blas::left(),
                      Common::Blas::upper(),
                      Common::Blas::no_trans(),
                      Common::Blas::non_unit(),
                      rows,
                      1,
                      1.,
                      A,
                      rows,
                      b,
                      1);
}

// Solves Ax = b, where AP = QR
void solve_qr_factorized(
    const double* QR, double* b, int rows, const int* permutations, const double* tau, double* work)
{
  // Calculate c = Q^T b;
  auto info = Common::Lapacke::dormqr(Common::Lapacke::row_major(), 'L', 'T', rows, 1, rows, QR, rows, tau, b, 1);
  if (info)
    DUNE_THROW(Dune::MathError, "Multiplication by Q^T failed");

  // Solve R x = c;
  solve_upper_triangular(QR, b, rows);

  // Undo permutations
  std::copy_n(b, rows, work);
  for (size_t ii = 0; ii < rows; ++ii)
    b[permutations[ii] - 1] = work[ii];
}


} // namespace LA
} // namespace GDT
} // namespace Dune
