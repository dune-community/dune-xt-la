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

#include <iostream>
#include <cmath>

#include <dune/common/exceptions.hh>

#include "solve_sym_tridiag_posdef.hh"

#include <dune/xt/common/lapacke.hh>

namespace Dune {
namespace XT {
namespace LA {


void qr(double* A, int rows, int cols, int* permutations, double* tau)
{
  auto info = Common::Lapacke::dgeqp3(Common::Lapacke::row_major(), rows, cols, A, cols, permutations, tau);
  if (info)
    DUNE_THROW(Dune::MathError, "QR factorization failed");
}


} // namespace LA
} // namespace GDT
} // namespace Dune
