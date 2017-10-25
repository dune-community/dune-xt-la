// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_XT_LA_CONTAINER_EYE_MATRIX_HH
#define DUNE_XT_LA_CONTAINER_EYE_MATRIX_HH

#include "pattern.hh"
#include "matrix-interface.hh"

namespace Dune {
namespace XT {
namespace LA {


template <class M>
typename std::enable_if<is_matrix<M>::value && M::sparse, M>::type eye_matrix(const size_t rows, const size_t cols)
{
  SparsityPatternDefault pattern(rows);
  for (size_t ii = 0; ii < std::min(rows, cols); ++ii)
    pattern.insert(ii, ii);
  // each row has to contain at least one non-zero entry
  for (size_t ii = std::min(rows, cols); ii < std::max(rows, cols); ++ii)
    pattern.insert(ii, 0);
  M mat(rows, cols, pattern);
  for (size_t ii = 0; ii < std::min(rows, cols); ++ii)
    mat.set_entry(ii, ii, 1);
  return mat;
}


template <class M>
typename std::enable_if<is_matrix<M>::value && !M::sparse, M>::type eye_matrix(const size_t rows, const size_t cols)
{
  M mat(rows, cols, 0.);
  for (size_t ii = 0; ii < std::min(rows, cols); ++ii)
    mat.set_entry(ii, ii, 1);
  return mat;
}


} // namespace LA
} // namespace XT
} // namespace Dune


#endif // DUNE_XT_LA_CONTAINER_EYE_MATRIX_HH