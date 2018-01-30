// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017 - 2018)

#ifndef DUNE_XT_LA_ALGORITHMS_SOLVE_TRIANGULAR_HH
#define DUNE_XT_LA_ALGORITHMS_SOLVE_TRIANGULAR_HH

#include <cstddef>

namespace Dune {
namespace XT {
namespace LA {

void solve_lower_triangular(const double* A, double* b, int rows);

void solve_lower_triangular_transposed(const double* A, double* b, int rows);

void solve_upper_triangular(const double* A, double* b, int rows);

/**  \brief Solves linear equation A*x = b for a matrix A that has been QR factorized by XT::LA::qr
  * \param[in] A triangular matrix
  * \param[in/out] b vector containing rhs, overwritten with solution vector x
  * \param[in] rows Number of rows of A
  * \param[in] cols Number of cols of A
  * \param[in] permutations Permutations from the QR factorization with column pivoting
  * \param[in] tau Multipliers from the QR factorization
  * \attention This function depends on LAPACKE. If LAPACKE is not found an error is thrown.
  */
double
solve_qr_factorized(const double* A, double* b, int rows, const int* permutations, const double* tau, double* work);


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_ALGORITHMS_SOLVE_TRIANGULAR_HH
