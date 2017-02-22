// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012 - 2014, 2016)
//   Rene Milk       (2013 - 2016)
//   Sven Kaulmann   (2014)

#ifndef DUNE_XT_LA_CONTAINER_PATTERN_HH
#define DUNE_XT_LA_CONTAINER_PATTERN_HH

#include <cstddef>
#include <vector>
#include <set>

namespace Dune {
namespace XT {
namespace LA {

class SparsityPatternDefault
{
private:
  typedef std::vector<std::vector<size_t>> BaseType;

public:
  typedef BaseType::value_type InnerType;
  typedef typename BaseType::const_iterator ConstOuterIteratorType;

  explicit SparsityPatternDefault(const size_t _size = 0);

  size_t size() const;

  InnerType& inner(const size_t ii);

  const InnerType& inner(const size_t ii) const;

  ConstOuterIteratorType begin() const;

  ConstOuterIteratorType end() const;

  bool operator==(const SparsityPatternDefault& other) const;

  bool operator!=(const SparsityPatternDefault& other) const;

  void insert(const size_t outer_index, const size_t inner_index);

  void sort(const size_t outer_index);

  void sort();

private:
  BaseType vector_of_vectors_;
}; // class SparsityPatternDefault

} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_CONTAINER_PATTERN_HH
