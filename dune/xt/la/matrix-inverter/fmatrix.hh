// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_XT_LA_MATRIX_INVERTER_FMATRIX_HH
#define DUNE_XT_LA_MATRIX_INVERTER_FMATRIX_HH

#include <dune/common/typetraits.hh>

#include <dune/xt/common/fmatrix.hh>
#include <dune/xt/la/exceptions.hh>
#include <dune/xt/la/matrix-inverter.hh>

#include "internal/base.hh"

namespace Dune {
namespace XT {
namespace LA {


template <class K, int ROWS, int COLS>
class MatrixInverterOptions<FieldMatrix<K, ROWS, COLS>>
{
public:
  static std::vector<std::string> types()
  {
    return {"direct"};
  }

  static Common::Configuration options(const std::string type = "")
  {
    const std::string actual_type = type.empty() ? types()[0] : type;
    internal::ensure_matrix_inverter_type(actual_type, types());
    Common::Configuration opts = internal::default_matrix_inverter_options();
    opts["type"] = actual_type;
    return opts;
  }
}; // class MatrixInverterOptions<EigenDenseMatrix<S>>


template <class K, int ROWS, int COLS>
class MatrixInverter<FieldMatrix<K, ROWS, COLS>> : public internal::MatrixInverterBase<FieldMatrix<K, ROWS, COLS>>
{
  using BaseType = internal::MatrixInverterBase<FieldMatrix<K, ROWS, COLS>>;

public:
  using MatrixType = typename BaseType::MatrixType;

  template <class... Args>
  explicit MatrixInverter(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
    const auto type = options_.template get<std::string>("type");
    const XT::Common::Configuration default_opts = MatrixInverterOptions<MatrixType>::options(type);
    if (!options_.get("delay_computation", default_opts.get<bool>("delay_computation")))
      compute();
  }

  void compute() override final
  {
    const auto type = options_.template get<std::string>("type");
    if (type == "direct") {
      inverse_ = std::make_unique<MatrixType>(matrix_);
      inverse_->invert();
    } else
      DUNE_THROW(Common::Exceptions::internal_error,
                 "Given type '" << type
                                << "' is none of MatrixInverterOptions<FieldMatrix<K, ROWS, COLS>>::types(), and "
                                   "internal::MatrixInverterBase promised to check this!"
                                << "\n\nThese are the available types:\n\n"
                                << MatrixInverterOptions<MatrixType>::types());

    this->post_checks();
  } // ... compute(...)

protected:
  using BaseType::matrix_;
  using BaseType::options_;
  using BaseType::inverse_;
}; // class MatrixInverter<EigenDenseMatrix<...>>


} // namespace Dune
} // namespace XT
} // namespace LA


#endif // DUNE_XT_LA_MATRIX_INVERTER_FMATRIX_HH