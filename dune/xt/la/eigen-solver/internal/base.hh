// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_XT_LA_EIGEN_SOLVER_INTERNAL_BASE_HH
#define DUNE_XT_LA_EIGEN_SOLVER_INTERNAL_BASE_HH

#include <algorithm>
#include <functional>
#include <memory>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/common/vector.hh>
#include <dune/xt/la/container/conversion.hh>
#include <dune/xt/la/container/matrix-interface.hh>
#include <dune/xt/la/exceptions.hh>
#include <dune/xt/la/matrix-inverter.hh>

namespace Dune {
namespace XT {
namespace LA {


// forward
template <class MatrixType>
class EigenSolverOptions;


namespace internal {


static void ensure_eigen_solver_type(const std::string& type, const std::vector<std::string>& available_types)
{
  bool contained = false;
  for (const auto& tp : available_types)
    if (type == tp)
      contained = true;
  if (!contained)
    DUNE_THROW(Exceptions::eigen_solver_failed_bc_it_was_not_set_up_correctly,
               "Given type '" << type << "' is not one of the available types: " << available_types);
} // ... ensure_type(...)


static Common::Configuration default_eigen_solver_options()
{
  Common::Configuration opts;
  opts["compute_eigenvalues"] = "true";
  opts["compute_eigenvectors"] = "true";
  opts["check_for_inf_nan"] = "true";
  opts["real_tolerance"] = "1e-15"; // is only used if the respective ensure_... is negative
  opts["ensure_real_eigenvalues"] = "-1"; // if positive, this is the check tolerance
  opts["ensure_positive_eigenvalues"] = "-1"; // if positive, this is the check tolerance
  opts["ensure_negative_eigenvalues"] = "-1"; // if positive, this is the check tolerance
  opts["ensure_real_eigenvectors"] = "-1"; // if positive, this is the check tolerance
  opts["check_eigendecomposition"] = "1e-10"; // disabled, if negative
  opts["check_real_eigendecomposition"] = "-1"; // if positive, this is the check tolerance
  return opts;
} // ... default_eigen_solver_options(...)


template <class MatrixImp, class FieldImp, class RealMatrixImp, class ComplexMatrixImp>
class EigenSolverBase
{
  static_assert(is_matrix<MatrixImp>::value || XT::Common::is_matrix<MatrixImp>::value, "");
  static_assert(is_matrix<RealMatrixImp>::value || XT::Common::is_matrix<RealMatrixImp>::value, "");
  static_assert(is_matrix<ComplexMatrixImp>::value || XT::Common::is_matrix<ComplexMatrixImp>::value, "");
  static_assert((is_matrix<MatrixImp>::value && is_matrix<RealMatrixImp>::value && is_matrix<ComplexMatrixImp>::value)
                    || (XT::Common::is_matrix<MatrixImp>::value && XT::Common::is_matrix<RealMatrixImp>::value
                        && XT::Common::is_matrix<ComplexMatrixImp>::value),
                "");
  using ThisType = EigenSolverBase<MatrixImp, FieldImp, RealMatrixImp, ComplexMatrixImp>;

public:
  using MatrixType = MatrixImp;
  using RealType = XT::Common::field_t<FieldImp>;
  using RealMatrixType = RealMatrixImp;
  using ComplexMatrixType = ComplexMatrixImp;

  EigenSolverBase(const MatrixType& matrix, const std::string& type = "")
    : matrix_(matrix)
    , options_(EigenSolverOptions<MatrixType>::options(type))
    , computed_(false)
    , eigenvalues_(nullptr)
    , real_eigenvalues_(nullptr)
    , eigenvectors_(nullptr)
    , real_eigenvectors_(nullptr)
  {
    pre_checks();
  }

  EigenSolverBase(const MatrixType& matrix, const Common::Configuration opts)
    : matrix_(matrix)
    , options_(opts)
    , computed_(false)
    , eigenvalues_(nullptr)
    , real_eigenvalues_(nullptr)
    , eigenvectors_(nullptr)
    , real_eigenvectors_(nullptr)
  {
    pre_checks();
  }

  virtual ~EigenSolverBase() = default;

  const Common::Configuration& options() const
  {
    return options_;
  }

  const MatrixType& matrix() const
  {
    return matrix_;
  }

protected:
  /**
   * \brief     Does the actual computation.
   * \attention The implementor has to fill the appropriate members!
   * \note      The implementor can assume that the given options_ contain a valid 'type' and all default keys.
   * \nte       The implementor does not need to guard against multiple calls of this method.
   */
  virtual void compute() const = 0;

public:
  const std::vector<XT::Common::complex_t<RealType>>& eigenvalues() const
  {
    compute_and_check();
    if (eigenvalues_)
      return *eigenvalues_;
    else if (options_.get<bool>("compute_eigenvalues"))
      DUNE_THROW(Common::Exceptions::internal_error, "The eigenvalues_ member is not filled after calling compute()!");
    else
      DUNE_THROW(Common::Exceptions::internal_error,
                 "Do not call eigenvalues() if 'compute_eigenvalues' is false!\n\nThese were the given options:\n\n"
                     << options_);
  } // ... eigenvalues(...)

  const std::vector<RealType>& real_eigenvalues() const
  {
    compute_and_check();
    if (eigenvalues_) {
      if (!real_eigenvalues_)
        compute_real_eigenvalues();
    } else if (options_.get<bool>("compute_eigenvalues"))
      DUNE_THROW(Common::Exceptions::internal_error, "The eigenvalues_ member is not filled after calling compute()!");
    else
      DUNE_THROW(
          Common::Exceptions::internal_error,
          "Do not call real_eigenvalues() if 'compute_eigenvalues' is false!\n\nThese were the given options:\n\n"
              << options_);
    assert(real_eigenvalues_ && "These have to exist after compute_real_eigenvalues()!");
    return *real_eigenvalues_;
  } // ... min_eigenvalues(...)

  std::vector<RealType> min_eigenvalues(const size_t num_evs = std::numeric_limits<size_t>::max()) const
  {
    compute_and_check();
    if (eigenvalues_) {
      if (!real_eigenvalues_)
        compute_real_eigenvalues();
    } else if (options_.get<bool>("compute_eigenvalues"))
      DUNE_THROW(Common::Exceptions::internal_error, "The eigenvalues_ member is not filled after calling compute()!");
    else
      DUNE_THROW(Common::Exceptions::internal_error,
                 "Do not call min_eigenvalues() if 'compute_eigenvalues' is false!\n\nThese were the given options:\n\n"
                     << options_);
    assert(real_eigenvalues_ && "These have to exist after compute_real_eigenvalues()!");
    std::vector<RealType> evs = *real_eigenvalues_;
    std::sort(evs.begin(), evs.end(), [](const RealType& a, const RealType& b) { return a < b; });
    evs.resize(std::min(evs.size(), num_evs));
    return evs;
  } // ... min_eigenvalues(...)

  std::vector<RealType> max_eigenvalues(const size_t num_evs = std::numeric_limits<size_t>::max()) const
  {
    compute_and_check();
    if (eigenvalues_) {
      if (!real_eigenvalues_)
        compute_real_eigenvalues();
    } else if (options_.get<bool>("compute_eigenvalues"))
      DUNE_THROW(Common::Exceptions::internal_error, "The eigenvalues_ member is not filled after calling compute()!");
    else
      DUNE_THROW(Common::Exceptions::internal_error,
                 "Do not call max_eigenvalues() if 'compute_eigenvalues' is false!\n\nThese were the given options:\n\n"
                     << options_);
    assert(real_eigenvalues_ && "These have to exist after compute_real_eigenvalues()!");
    std::vector<RealType> evs = *real_eigenvalues_;
    std::sort(evs.begin(), evs.end(), [](const RealType& a, const RealType& b) { return a > b; });
    evs.resize(std::min(evs.size(), num_evs));
    return evs;
  } // ... max_eigenvalues(...)

  const ComplexMatrixType& eigenvectors() const
  {
    compute_and_check();
    if (eigenvectors_)
      return *eigenvectors_;
    else if (options_.get<bool>("compute_eigenvectors"))
      DUNE_THROW(Common::Exceptions::internal_error, "The eigenvectors_ member is not filled after calling compute()!");
    else
      DUNE_THROW(Common::Exceptions::internal_error,
                 "Do not call eigenvectors() if 'compute_eigenvectors' is false!\n\nThese were the given options:\n\n"
                     << options_);
  } // ... eigenvectors(...)

  const RealMatrixType& real_eigenvectors() const
  {
    compute_and_check();
    if (eigenvectors_) {
      if (!real_eigenvectors_)
        compute_real_eigenvectors();
    } else if (options_.get<bool>("compute_eigenvectors"))
      DUNE_THROW(Common::Exceptions::internal_error, "The eigenvectors_ member is not filled after calling compute()!");
    else
      DUNE_THROW(
          Common::Exceptions::internal_error,
          "Do not call real_eigenvectors() if 'compute_eigenvectors' is false!\n\nThese were the given options:\n\n"
              << options_);
    assert(real_eigenvectors_ && "These have to exist after compute_real_eigenvectors()!");
    return *real_eigenvectors_;
  } // ... real_eigenvectors(...)

protected:
  void compute_and_check() const
  {
    if (!computed_) {
      compute();
      post_checks();
    }
    computed_ = true;
  }

  void pre_checks()
  {
    // check options
    if (!options_.has_key("type"))
      DUNE_THROW(Exceptions::eigen_solver_failed_bc_it_was_not_set_up_correctly,
                 "Missing 'type' in given options!"
                     << "\n\nThese were the given options:\n\n"
                     << options_);
    internal::ensure_eigen_solver_type(options_.get<std::string>("type"), EigenSolverOptions<MatrixType>::types());
    const Common::Configuration default_opts =
        EigenSolverOptions<MatrixType>::options(options_.get<std::string>("type"));
    for (const std::string& default_key : default_opts.getValueKeys()) {
      if (!options_.has_key(default_key))
        options_[default_key] = default_opts.get<std::string>(default_key);
    }
    if (options_.get<double>("ensure_positive_eigenvalues") > 0
        && options_.get<double>("ensure_negative_eigenvalues") > 0)
      DUNE_THROW(Exceptions::eigen_solver_failed_bc_it_was_not_set_up_correctly,
                 "It does not make sense to ensure positive and negative eigenvalues!");
    if (!options_.get<bool>("compute_eigenvalues") && (options_.get<double>("ensure_real_eigenvalues") > 0
                                                       || options_.get<double>("ensure_positive_eigenvalues") > 0
                                                       || options_.get<double>("ensure_negative_eigenvalues") > 0))
      options_["compute_eigenvalues"] = true;
    if (options_.get<double>("check_real_eigendecomposition") > 0) {
      if (!options_.get<bool>("compute_eigenvalues"))
        options_["compute_eigenvalues"] = "true";
      if (options_.get<double>("ensure_real_eigenvalues") <= 0)
        options_["ensure_real_eigenvalues"] = options_["check_real_eigendecomposition"];
    }
    // check matrix
    check_size(matrix_);
    if (options_.get<bool>("check_for_inf_nan") && contains_inf_or_nan(matrix_)) {
      DUNE_THROW(Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements,
                 "Given matrix contains inf or nan and you requested checking. To disable this check set "
                 "'check_for_inf_nan' to false in the options."
                     << "\n\nThese were the given options:\n\n"
                     << options_
                     << "\nThis was the given matrix:\n\n"
                     << matrix_);
    }
  } // ... pre_checks(...)

  void post_checks() const
  {
    if (options_.get<bool>("compute_eigenvalues") && !eigenvalues_)
      DUNE_THROW(Common::Exceptions::internal_error, "The eigenvalues_ member is not filled after calling compute()!");
    if (options_.get<bool>("compute_eigenvectors") && !eigenvectors_)
      DUNE_THROW(Common::Exceptions::internal_error, "The eigenvectors_ member is not filled after calling compute()!");
    if (options_.get<bool>("check_for_inf_nan")) {
      if (eigenvalues_ && contains_inf_or_nan(*eigenvalues_))
        DUNE_THROW(Exceptions::eigen_solver_failed_bc_result_contained_inf_or_nan,
                   "Computed eigenvalues contain inf or nan and you requested checking. To disable this check set "
                   "'check_for_inf_nan' to false in the options."
                       << "\n\nThese were the given options:\n\n"
                       << options_
                       << "\nThese are the computed eigenvalues:\n\n"
                       << *eigenvalues_);
      if (eigenvectors_ && contains_inf_or_nan(*eigenvectors_))
        DUNE_THROW(Exceptions::eigen_solver_failed_bc_result_contained_inf_or_nan,
                   "Computed eigenvectors contain inf or nan and you requested checking. To disable this check set "
                   "'check_for_inf_nan' to false in the options."
                       << "\n\nThese were the given options:\n\n"
                       << options_
                       << "\nThese are the computed eigenvectors:\n\n"
                       << *eigenvectors_);
    }
    const double ensure_real_eigenvalues = options_.get<double>("ensure_real_eigenvalues");
    const double ensure_positive_eigenvalues = options_.get<double>("ensure_positive_eigenvalues");
    const double ensure_negative_eigenvalues = options_.get<double>("ensure_negative_eigenvalues");
    const double check_real_eigendecomposition = options_.get<double>("check_real_eigendecomposition");
    if (ensure_real_eigenvalues > 0 || ensure_positive_eigenvalues > 0 || ensure_negative_eigenvalues > 0
        || check_real_eigendecomposition > 0)
      compute_real_eigenvalues();
    if (ensure_positive_eigenvalues > 0) {
      assert(real_eigenvalues_ && "This must not happen after compute_real_eigenvalues()!");
      for (const auto& ev : *real_eigenvalues_) {
        if (ev < ensure_positive_eigenvalues)
          DUNE_THROW(Exceptions::eigen_solver_failed_bc_eigenvalues_are_not_positive_as_requested,
                     "These were the given options:\n\n"
                         << options_
                         << "\nThese are the computed eigenvectors:\n\n"
                         << *eigenvectors_);
      }
    }
    if (ensure_negative_eigenvalues > 0) {
      assert(real_eigenvalues_ && "This must not happen after compute_real_eigenvalues()!");
      for (const auto& ev : *real_eigenvalues_) {
        if (ev > -1 * ensure_negative_eigenvalues)
          DUNE_THROW(Exceptions::eigen_solver_failed_bc_eigenvalues_are_not_negative_as_requested,
                     "These were the given options:\n\n"
                         << options_
                         << "\nThese are the computed eigenvectors:\n\n"
                         << *eigenvectors_);
      }
    }
    if (options_.get<double>("ensure_real_eigenvectors") > 0 || check_real_eigendecomposition > 0)
      compute_real_eigenvectors();
    const double check_eigendecomposition = options_.get<double>("check_eigendecomposition");
    if (check_eigendecomposition > 0)
      complex_eigendecomposition_helper<>::check(*this, check_eigendecomposition);
    if (check_real_eigendecomposition > 0)
      assert_eigendecomposition(matrix_, *real_eigenvalues_, *real_eigenvectors_, check_real_eigendecomposition);
  } // ... post_checks(...)

  void compute_real_eigenvalues() const
  {
    assert(eigenvalues_ && "This should not happen!");
    if (!real_eigenvalues_) {
      const double ensure_real_eigenvalues = options_.get<double>("ensure_real_eigenvalues");
      const double tolerance =
          (ensure_real_eigenvalues > 0) ? ensure_real_eigenvalues : options_.get<double>("real_tolerance");
      real_eigenvalues_ = std::make_unique<std::vector<RealType>>(eigenvalues_->size());
      for (size_t ii = 0; ii < eigenvalues_->size(); ++ii) {
        const auto& complex_ev = (*eigenvalues_)[ii];
        if (std::abs(complex_ev.imag()) > tolerance)
          DUNE_THROW(Exceptions::eigen_solver_failed_bc_eigenvalues_are_not_real_as_requested,
                     "These were the given options:\n\n"
                         << options_
                         << "\nThese are the computed eigenvalues:\n\n"
                         << *eigenvalues_);
        (*real_eigenvalues_)[ii] = complex_ev.real();
      }
    }
  } // ... compute_real_eigenvalues(...)

  template <bool is_common_matrix = XT::Common::is_matrix<MatrixType>::value, class T = MatrixType>
  struct real_eigenvectors_helper
  {
  };

  template <class T>
  struct real_eigenvectors_helper<true, T>
  {
    static void compute(const ThisType& self, const double& tolerance)
    {
      using RM = XT::Common::MatrixAbstraction<RealMatrixType>;
      using CM = XT::Common::MatrixAbstraction<ComplexMatrixType>;
      const size_t rows = CM::rows(*self.eigenvectors_);
      const size_t cols = CM::cols(*self.eigenvectors_);
      self.real_eigenvectors_ = std::make_unique<RealMatrixType>(RM::create(rows, cols));
      for (size_t ii = 0; ii < rows; ++ii)
        for (size_t jj = 0; jj < cols; ++jj) {
          const auto complex_value = CM::get_entry(*self.eigenvectors_, ii, jj);
          if (std::abs(complex_value.imag()) > tolerance)
            DUNE_THROW(Exceptions::eigen_solver_failed_bc_eigenvectors_are_not_real_as_requested,
                       "These were the given options:\n\n"
                           << self.options_
                           << "\nThese are the computed eigenvectors:\n\n"
                           << *self.eigenvectors_);
          RM::set_entry(*self.real_eigenvectors_, ii, jj, complex_value.real());
        }
    }
  }; // real_eigenvectors_helper<true, ...>

  template <class T>
  struct real_eigenvectors_helper<false, T>
  {
    static void compute(ThisType& self, const double& tolerance)
    {
      if (RealMatrixType::sparse) {
        const size_t rows = self.eigenvectors_->rows();
        const size_t cols = self.eigenvectors_->cols();
        const auto pattern = self.eigenvectors_->pattern();
        self.real_eigenvectors_ = std::make_unique<RealMatrixType>(rows, cols, pattern);
        for (size_t ii = 0; ii < rows; ++ii)
          for (size_t jj : pattern.inner(ii)) {
            const auto complex_value = self.eigenvectors_->get_entry(ii, jj);
            if (std::abs(complex_value.imag()) > tolerance)
              DUNE_THROW(Exceptions::eigen_solver_failed_bc_eigenvectors_are_not_real_as_requested,
                         "These were the given options:\n\n"
                             << self.options_
                             << "\nThese are the computed eigenvectors:\n\n"
                             << *self.eigenvectors_);
            self.real_eigenvectors_->set_entry(ii, jj, complex_value.real());
          }
      } else {
        const size_t rows = self.eigenvectors_->rows();
        const size_t cols = self.eigenvectors_->cols();
        self.real_eigenvectors_ = std::make_unique<RealMatrixType>(rows, cols);
        for (size_t ii = 0; ii < rows; ++ii)
          for (size_t jj = 0; jj < cols; ++jj) {
            const auto complex_value = self.eigenvectors_->get_entry(ii, jj);
            if (std::abs(complex_value.imag()) > tolerance)
              DUNE_THROW(Exceptions::eigen_solver_failed_bc_eigenvectors_are_not_real_as_requested,
                         "These were the given options:\n\n"
                             << self.options_
                             << "\nThese are the computed eigenvectors:\n\n"
                             << *self.eigenvectors_);
            self.real_eigenvectors_->set_entry(ii, jj, complex_value.real());
          }
      }
    }
  }; // real_eigenvectors_helper<false, ...>

  void compute_real_eigenvectors() const
  {
    assert(eigenvectors_ && "This should not happen!");
    if (!real_eigenvectors_) {
      const double ensure_real_eigenvectors = options_.get<double>("ensure_real_eigenvectors");
      const double tolerance =
          (ensure_real_eigenvectors > 0) ? ensure_real_eigenvectors : options_.get<double>("real_tolerance");
      real_eigenvectors_helper<>::compute(*this, tolerance);
    }
  }

  template <class A, class B, class C>
  void
  assert_eigendecomposition(const A& mat, const B& eigenvalues, const C& eigenvectors, const double& tolerance) const
  {
    const size_t rows = Common::get_matrix_rows(mat);
    const size_t cols = Common::get_matrix_cols(mat);
    auto eigenvectors_inv = Common::create<C>(cols, rows);
    try {
      eigenvectors_inv = invert_matrix(eigenvectors);
    } catch (const Exceptions::matrix_invert_failed& ee) {
      DUNE_THROW(Exceptions::eigen_solver_failed,
                 "The computed matrix of eigenvectors is not invertible!"
                     << "\n\nmatrix = "
                     << matrix_
                     << "\n\noptions: "
                     << options_
                     << "\n\neigenvalues = "
                     << eigenvalues
                     << "\n\neigenvectors = "
                     << eigenvectors
                     << "\n\nThis was the original error: "
                     << ee.what());
    } catch (const Dune::Exception& ee) {
      DUNE_THROW(Exceptions::eigen_solver_failed,
                 "An error occured during inversion of the matrix of eigenvectors!"
                     << "\n\nmatrix = "
                     << matrix_
                     << "\n\noptions: "
                     << options_
                     << "\n\neigenvalues = "
                     << eigenvalues
                     << "\n\neigenvectors = "
                     << eigenvectors
                     << "\n\nThis was the original error: "
                     << ee.what());
    } catch (...) {
      DUNE_THROW(Exceptions::eigen_solver_failed,
                 "An unknown error occured during inversion of the matrix of eigenvectors!"
                     << "\n\nmatrix = "
                     << matrix_
                     << "\n\noptions: "
                     << options_
                     << "\n\neigenvalues = "
                     << eigenvalues
                     << "\n\neigenvectors = "
                     << eigenvectors);
    }
    auto eigenvalue_matrix = Common::MatrixAbstraction<C>::create(rows, cols, 0.);
    for (size_t ii = 0; ii < rows; ++ii)
      Common::set_matrix_entry(eigenvalue_matrix, ii, ii, eigenvalues[ii]);
    const auto decomposition_error = (eigenvectors * (eigenvalue_matrix * eigenvectors_inv)) - mat;
    for (size_t ii = 0; ii < rows; ++ii)
      for (size_t jj = 0; jj < cols; ++jj)
        if (std::abs(Common::get_matrix_entry(decomposition_error, ii, jj)) > tolerance)
          DUNE_THROW(Exceptions::eigen_solver_failed_bc_result_is_not_an_eigendecomposition,
                     "\n\nmatrix = " << matrix_ << "\n\noptions: " << options_ << "\n\neigenvalues (lambda)= "
                                     << eigenvalues
                                     << "\n\neigenvectors (T) = "
                                     << eigenvectors
                                     << "\n\n(T * (lambda * T^-1)) - matrix = "
                                     << (eigenvectors * (eigenvalue_matrix * eigenvectors_inv)) - mat);
  } // ... assert_eigendecomposition(...)

  template <bool upcast_required = !std::is_same<MatrixType, ComplexMatrixType>::value, bool anything = true>
  struct complex_eigendecomposition_helper;

  template <bool anything>
  struct complex_eigendecomposition_helper<true, anything>
  {
    static void check(const ThisType& self, const double& tolerance)
    {
      self.assert_eigendecomposition(
          convert_to<ComplexMatrixType>(self.matrix_), *self.eigenvalues_, *self.eigenvectors_, tolerance);
    }
  };

  template <bool anything>
  struct complex_eigendecomposition_helper<false, anything>
  {
    static void check(const ThisType& self, const double& tolerance)
    {
      self.assert_eigendecomposition(self.matrix_, *self.eigenvalues_, *self.eigenvectors_, tolerance);
    }
  };

  template <class M>
  void check_size(const MatrixInterface<M>& mat) const
  {
    if (mat.rows() != mat.cols())
      DUNE_THROW(Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements,
                 "Matrix has to be square, is " << mat.rows() << "x" << mat.cols() << "!");
  }

  template <class M>
  typename std::enable_if<XT::Common::is_matrix<M>::value && !is_matrix<M>::value, void>::type
  check_size(const M& mat) const
  {
    using Mat = XT::Common::MatrixAbstraction<M>;
    if (Mat::rows(mat) != Mat::cols(mat))
      DUNE_THROW(Exceptions::eigen_solver_failed_bc_data_did_not_fulfill_requirements,
                 "Matrix has to be square, is " << Mat::rows(mat) << "x" << Mat::cols(mat) << "!");
  }

  template <class T>
  bool contains_inf_or_nan(const std::vector<T>& vec) const
  {
    for (const auto& element : vec)
      if (XT::Common::isinf(element) || XT::Common::isnan(element))
        return true;
    return false;
  }

  template <class M>
  bool contains_inf_or_nan(const MatrixInterface<M>& mat) const
  {
    return !mat.valid();
  }

  template <class M>
  typename std::enable_if<XT::Common::is_matrix<M>::value && !is_matrix<M>::value, bool>::type
  contains_inf_or_nan(const M& mat) const
  {
    using Mat = XT::Common::MatrixAbstraction<M>;
    for (size_t ii = 0; ii < Mat::rows(mat); ++ii)
      for (size_t jj = 0; jj < Mat::cols(mat); ++jj) {
        const auto value = Mat::get_entry(mat, ii, jj);
        if (XT::Common::isinf(value) || XT::Common::isnan(value))
          return true;
      }
    return false;
  } // ... contains_inf_or_nan(...)

  const MatrixType& matrix_;
  mutable Common::Configuration options_;
  mutable bool computed_;
  mutable std::unique_ptr<std::vector<XT::Common::complex_t<RealType>>> eigenvalues_;
  mutable std::unique_ptr<std::vector<RealType>> real_eigenvalues_;
  mutable std::unique_ptr<ComplexMatrixType> eigenvectors_;
  mutable std::unique_ptr<RealMatrixType> real_eigenvectors_;
}; // class EigenSolverBase


} // namespace internal
} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_EIGEN_SOLVER_INTERNAL_BASE_HH