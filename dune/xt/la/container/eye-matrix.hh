// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)

#ifndef DUNE_XT_LA_CONTAINER_EYE_MATRIX_HH
#define DUNE_XT_LA_CONTAINER_EYE_MATRIX_HH

#include "pattern.hh"
#include "matrix-interface.hh"

namespace Dune {
namespace XT {
namespace LA {
namespace internal {


SparsityPatternDefault create_eye_pattern(const size_t rows, const size_t cols)
{
  SparsityPatternDefault pattern(rows);
  for (size_t ii = 0; ii < std::min(rows, cols); ++ii)
    pattern.insert(ii, ii);
  // each row has to contain at least one non-zero entry
  for (size_t ii = std::min(rows, cols); ii < std::max(rows, cols); ++ii)
    pattern.insert(ii, 0);
  return pattern;
}

template <class M>
typename std::enable_if<is_matrix<M>::value, void>::type set_diagonal_to_one(M& mat)
{
  for (size_t ii = 0; ii < std::min(mat.rows(), mat.cols()); ++ii)
    mat.set_entry(ii, ii, 1.);
}


} // namespace internal


template <class M>
typename std::enable_if<is_matrix<M>::value, M>::type eye_matrix(const size_t rows, const size_t cols)
{
  M mat = M::sparse ? M(rows, cols, internal::create_eye_pattern(rows, cols)) : M(rows, cols, 0.);
  internal::set_diagonal_to_one(mat);
  return mat;
}

template <class M>
typename std::enable_if<is_matrix<M>::value, std::unique_ptr<M>>::type eye_matrix_ptr(const size_t rows,
                                                                                      const size_t cols)
{
  auto mat = M::sparse ? std::make_unique<M>(rows, cols, internal::create_eye_pattern(rows, cols))
                       : std::make_unique<M>(rows, cols, 0.);
  internal::set_diagonal_to_one(*mat);
  return mat;
}

template <class M>
typename std::enable_if<Common::is_matrix<M>::value && !is_matrix<M>::value, M>::type eye_matrix(const size_t rows,
                                                                                                 const size_t cols)
{
  using Abstraction = Common::MatrixAbstraction<M>;
  auto mat = Abstraction::create(rows, cols, 0.);
  for (size_t ii = 0; ii < std::min(rows, cols); ++ii)
    Abstraction::set_entry(mat, ii, ii, 1);
  return mat;
}

template <class M>
typename std::enable_if<Common::is_matrix<M>::value && !is_matrix<M>::value, std::unique_ptr<M>>::type
eye_matrix_ptr(const size_t rows, const size_t cols)
{
  using Abstraction = Common::MatrixAbstraction<M>;
  auto mat = Abstraction::create_dynamic(rows, cols, 0.);
  for (size_t ii = 0; ii < std::min(rows, cols); ++ii)
    Abstraction::set_entry(*mat, ii, ii, 1);
  return mat;
}

template <class M>
typename std::enable_if<Common::is_matrix<M>::value, M>::type eye_matrix(const size_t size)
{
  return eye_matrix<M>(size, size);
}

template <class M>
typename std::enable_if<Common::is_matrix<M>::value, std::unique_ptr<M>>::type eye_matrix_ptr(const size_t size)
{
  return eye_matrix_ptr<M>(size, size);
}


} // namespace LA
} // namespace XT
} // namespace Dune


#endif // DUNE_XT_LA_CONTAINER_EYE_MATRIX_HH
