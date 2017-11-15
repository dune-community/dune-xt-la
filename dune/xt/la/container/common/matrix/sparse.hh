// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Barbara Verfürth (2015)
//   Felix Schindler  (2013 - 2017)
//   Rene Milk        (2014 - 2016)
//   Tobias Leibner   (2014, 2016 - 2017)

#ifndef DUNE_XT_LA_CONTAINER_COMMON_MATRIX_SPARSE_HH
#define DUNE_XT_LA_CONTAINER_COMMON_MATRIX_SPARSE_HH

#include <dune/xt/la/container/interfaces.hh>
#include <dune/xt/la/container/pattern.hh>

#include "../vector/sparse.hh"
#include "dense.hh"

namespace Dune {
namespace XT {
namespace LA {


enum class SparseFormat
{
  csr,
  csc
};

// forwards
template <class ScalarImp, SparseFormat sparse_format>
class CommonSparseMatrix;

// forwards
template <class DenseMatrixImp, class SparseMatrixImp>
class CommonSparseOrDenseMatrix;


namespace internal {


template <class ScalarImp, SparseFormat sparse_format>
class CommonSparseMatrixTraits
{
public:
  typedef typename Dune::FieldTraits<ScalarImp>::field_type ScalarType;
  typedef typename Dune::FieldTraits<ScalarImp>::real_type RealType;
  static const Backends backend_type = Backends::common_dense;
  static const Backends vector_type = Backends::common_dense;
  typedef std::vector<ScalarImp> EntriesVectorType;
  typedef std::vector<size_t> IndexVectorType;
  typedef typename Common::FloatCmp::DefaultEpsilon<ScalarType>::Type EpsType;
  typedef CommonSparseMatrix<ScalarImp, sparse_format> derived_type;
  static const constexpr bool sparse = true;
};

template <class DenseMatrixImp, class SparseMatrixImp>
class CommonSparseOrDenseMatrixTraits
{
public:
  typedef DenseMatrixImp DenseMatrixType;
  typedef SparseMatrixImp SparseMatrixType;
  typedef typename DenseMatrixType::ScalarType ScalarType;
  typedef typename DenseMatrixType::RealType RealType;
  static const Backends backend_type = Backends::common_dense;
  static const Backends vector_type = Backends::common_dense;
  typedef typename Common::FloatCmp::DefaultEpsilon<ScalarType>::Type EpsType;
  typedef CommonSparseOrDenseMatrix<DenseMatrixType, SparseMatrixType> derived_type;
  static const constexpr bool sparse = false;
};


} // namespace internal


/**
 * \brief A sparse matrix implementation of the MatrixInterface with row major memory layout.
 */
template <class ScalarImp = double, SparseFormat sparse_format = SparseFormat::csr>
class CommonSparseMatrix
    : public MatrixInterface<internal::CommonSparseMatrixTraits<ScalarImp, sparse_format>, ScalarImp>
{
  typedef CommonSparseMatrix<ScalarImp, sparse_format> ThisType;
  typedef MatrixInterface<internal::CommonSparseMatrixTraits<ScalarImp, sparse_format>, ScalarImp> MatrixInterfaceType;

public:
  typedef internal::CommonSparseMatrixTraits<ScalarImp, sparse_format> Traits;
  typedef typename Traits::EntriesVectorType EntriesVectorType;
  typedef typename Traits::IndexVectorType IndexVectorType;
  typedef typename Traits::ScalarType ScalarType;
  typedef typename Traits::RealType RealType;
  typedef typename Traits::EpsType EpsType;

  /**
  * \brief This is the constructor of interest which creates a sparse matrix.
  */
  CommonSparseMatrix(const size_t rr,
                     const size_t cc,
                     const SparsityPatternDefault& patt,
                     const size_t num_mutexes = 1,
                     const EpsType eps = Common::FloatCmp::DefaultEpsilon<ScalarType>::value() / 1000.)
    : num_rows_(rr)
    , num_cols_(cc)
    , entries_(std::make_shared<EntriesVectorType>())
    , row_pointers_(std::make_shared<IndexVectorType>(num_rows_ + 1, 0))
    , column_indices_(std::make_shared<IndexVectorType>())
    , mutexes_(num_mutexes > 0 ? std::make_shared<std::vector<std::mutex>>(num_mutexes) : nullptr)
    , eps_(eps)
    , unshareable_(false)
  {
    if (num_rows_ > 0 && num_cols_ > 0) {
      if (size_t(patt.size()) != num_rows_)
        DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                   "The size of the pattern (" << patt.size() << ") does not match the number of rows of this ("
                                               << num_rows_
                                               << ")!");
      for (size_t row = 0; row < num_rows_; ++row) {
        const auto& columns = patt.inner(row);
        const auto num_nonzero_entries_in_row = columns.size();
        assert(columns.size() <= num_cols_);
        row_pointers_->operator[](row + 1) = row_pointers_->operator[](row) + num_nonzero_entries_in_row;
        for (size_t kk = 0; kk < num_nonzero_entries_in_row; ++kk) {
#ifndef NDEBUG
          if (columns[kk] >= num_cols_)
            DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                       "The size of row " << row << " of the pattern does not match the number of columns of this ("
                                          << num_cols_
                                          << ")!");
#endif // NDEBUG
          column_indices_->push_back(columns[kk]);
        } // kk
        entries_->resize(column_indices_->size());
      } // row
    }
  } // CommonSparseMatrix(rr, cc, patt, num_mutexes)

  CommonSparseMatrix(const size_t rr = 0,
                     const size_t cc = 0,
                     const ScalarType& value = ScalarType(0),
                     const size_t num_mutexes = 1,
                     const EpsType eps = Common::FloatCmp::DefaultEpsilon<ScalarType>::value() / 1000.)
    : num_rows_(rr)
    , num_cols_(cc)
    , entries_(std::make_shared<EntriesVectorType>(
          XT::Common::FloatCmp::eq(value, ScalarType(0.), 0., eps / num_rows_) ? 0 : num_rows_ * num_cols_, value))
    , row_pointers_(std::make_shared<IndexVectorType>(num_rows_ + 1, 0))
    , column_indices_(std::make_shared<IndexVectorType>())
    , mutexes_(num_mutexes > 0 ? std::make_shared<std::vector<std::mutex>>(num_mutexes) : nullptr)
    , eps_(eps)
    , unshareable_(false)
  {
    if (XT::Common::FloatCmp::ne(value, ScalarType(0), 0., eps / num_rows_)) {
      IndexVectorType row_column_indices(num_cols_);
      for (size_t col = 0; col < num_cols_; ++col)
        row_column_indices[col] = col;
      for (size_t row = 0; row < num_rows_; ++row) {
        row_pointers_->operator[](row + 1) = (row + 1) * num_cols_;
        column_indices_->insert(column_indices_->end(), row_column_indices.begin(), row_column_indices.end());
      }
    }
  }

  CommonSparseMatrix(const ThisType& other)
    : num_rows_(other.num_rows_)
    , num_cols_(other.num_cols_)
    , entries_(other.unshareable_ ? std::make_shared<EntriesVectorType>(*other.entries_) : other.entries_)
    , row_pointers_(other.row_pointers_)
    , column_indices_(other.column_indices_)
    , mutexes_(other.unshareable_ ? std::make_shared<std::vector<std::mutex>>(other.mutexes_->size()) : other.mutexes_)
    , eps_(other.eps_)
    , unshareable_(false)
  {
  }

  template <class OtherMatrixType>
  explicit CommonSparseMatrix(
      const OtherMatrixType& mat,
      const typename std::enable_if<Common::MatrixAbstraction<OtherMatrixType>::is_matrix, bool>::type prune = false,
      const EpsType eps = Common::FloatCmp::DefaultEpsilon<ScalarType>::value() / 1000.,
      const size_t num_mutexes = 1)
    : num_rows_(Common::MatrixAbstraction<OtherMatrixType>::rows(mat))
    , num_cols_(Common::MatrixAbstraction<OtherMatrixType>::cols(mat))
    , entries_(std::make_shared<EntriesVectorType>())
    , row_pointers_(std::make_shared<IndexVectorType>(num_rows_ + 1))
    , column_indices_(std::make_shared<IndexVectorType>())
    , mutexes_(num_mutexes > 0 ? std::make_shared<std::vector<std::mutex>>(num_mutexes) : nullptr)
    , eps_(eps)
    , unshareable_(false)
  {
    size_t index = 0;
    for (size_t rr = 0; rr < num_rows_; ++rr) {
      (*row_pointers_)[rr] = index;
      for (size_t cc = 0; cc < num_cols_; ++cc) {
        const auto& value = Common::MatrixAbstraction<OtherMatrixType>::get_entry(mat, rr, cc);
        if (!prune || XT::Common::FloatCmp::ne(
                          value, ScalarType(0), 0., eps / Common::MatrixAbstraction<OtherMatrixType>::cols(mat))) {
          entries_->push_back(value);
          column_indices_->push_back(cc);
          ++index;
        }
      } // cc
    } // rr
    (*row_pointers_)[num_rows_] = column_indices_->size();
  } // CommonSparseMatrix(...)

  template <class OtherMatrixType>
  explicit CommonSparseMatrix(
      const OtherMatrixType& mat,
      const typename std::enable_if<Common::MatrixAbstraction<OtherMatrixType>::is_matrix, bool>::type prune,
      const size_t num_mutexes)
    : CommonSparseMatrix(mat,
                         prune,
                         Common::FloatCmp::DefaultEpsilon<ScalarType>::value()
                             / (1000. * Common::MatrixAbstraction<OtherMatrixType>::cols(mat)),
                         num_mutexes)
  {
  } // CommonSparseMatrix(...)

  template <class DenseMatrixImp>
  void copy_to_densematrix(DenseMatrixImp& ret) const
  {
    for (size_t rr = 0; rr < num_rows_; ++rr)
      for (size_t kk = row_pointers_->operator[](rr); kk < row_pointers_->operator[](rr + 1); ++kk)
        ret[rr][column_indices_->operator[](kk)] = entries_->operator[](kk);
  }

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      num_rows_ = other.num_rows_;
      num_cols_ = other.num_cols_;
      entries_ = other.unshareable_ ? std::make_shared<EntriesVectorType>(*other.entries_) : other.entries_;
      row_pointers_ = other.row_pointers_;
      column_indices_ = other.column_indices_;
      mutexes_ = other.mutexes_
                     ? (other.unshareable_ ? std::make_shared<std::vector<std::mutex>>(other.mutexes_->size())
                                           : other.mutexes_)
                     : nullptr;
      unshareable_ = false;
    }
    return *this;
  }

  template <class OtherMatrixImp>
  typename std::enable_if_t<XT::Common::MatrixAbstraction<OtherMatrixImp>::is_matrix, ThisType>&
  operator=(const OtherMatrixImp& other)
  {
    clear();
    typedef XT::Common::MatrixAbstraction<OtherMatrixImp> MatAbstrType;
    num_rows_ = MatAbstrType::rows(other);
    num_cols_ = MatAbstrType::cols(other);
    for (size_t rr = 0; rr < num_rows_; ++rr) {
      for (size_t cc = 0; cc < num_cols_; ++cc) {
        if (XT::Common::FloatCmp::ne(MatAbstrType::get_entry(other, rr, cc), 0., 0., eps_ / num_cols_)) {
          entries_->push_back(MatAbstrType::get_entry(other, rr, cc));
          column_indices_->push_back(cc);
        }
      } // cc
      (*row_pointers_)[rr + 1] = column_indices_->size();
    } // rr
    return *this;
  }

  void deep_copy(const ThisType& other)
  {
    ensure_uniqueness();
    num_rows_ = other.num_rows_;
    num_cols_ = other.num_cols_;
    *entries_ = *other.entries_;
    *row_pointers_ = *other.row_pointers_;
    *column_indices_ = *other.column_indices_;
  }

  void clear()
  {
    ensure_uniqueness();
    entries_->clear();
    std::fill(row_pointers_->begin(), row_pointers_->end(), 0);
    column_indices_->clear();
  }

  /// \name Required by ContainerInterface.
  /// \{
  inline ThisType copy() const
  {
    ThisType ret(*this);
    ret.entries_ = std::make_shared<EntriesVectorType>(*entries_);
    ret.row_pointers_ = std::make_shared<IndexVectorType>(*row_pointers_);
    ret.column_indices_ = std::make_shared<IndexVectorType>(*column_indices_);
    ret.mutexes_ = mutexes_ ? std::make_shared<std::vector<std::mutex>>(mutexes_->size()) : nullptr;
    return ret;
  }

  inline void scal(const ScalarType& alpha)
  {
    ensure_uniqueness();
    const internal::VectorLockGuard DUNE_UNUSED(guard)(mutexes_);
    std::transform(
        entries_->begin(), entries_->end(), entries_->begin(), std::bind1st(std::multiplies<ScalarType>(), alpha));
  }

  inline void axpy(const ScalarType& alpha, const ThisType& xx)
  {
    ensure_uniqueness();
    const internal::VectorLockGuard DUNE_UNUSED(guard)(mutexes_);
    assert(has_equal_shape(xx));
    const auto& xx_entries = *xx.entries_;
    for (size_t ii = 0; ii < entries_->size(); ++ii)
      entries_->operator[](ii) += alpha * xx_entries[ii];
  }

  inline bool has_equal_shape(const ThisType& other) const
  {
    return (rows() == other.rows()) && (cols() == other.cols());
  }

  /// \}
  /// \name Required by MatrixInterface.
  /// \{

  inline size_t rows() const
  {
    return num_rows_;
  }

  inline size_t cols() const
  {
    return num_cols_;
  }

  //! Matrix-Vector multiplication for arbitrary vectors that support operator[]
  template <class XX, class YY>
  inline std::enable_if_t<XT::Common::VectorAbstraction<XX>::is_vector && XT::Common::VectorAbstraction<YY>::is_vector,
                          void>
  mv(const XX& xx, YY& yy) const
  {
    std::fill(yy.begin(), yy.end(), ScalarType(0));
    const auto& entries = *entries_;
    const auto& row_pointers = *row_pointers_;
    const auto& column_indices = *column_indices_;
    for (size_t rr = 0; rr < num_rows_; ++rr) {
      auto& yy_rr = yy[rr];
      const size_t end = row_pointers[rr + 1];
      for (size_t kk = row_pointers[rr]; kk < end; ++kk)
        yy_rr += entries[kk] * xx[column_indices[kk]];
    }
  }

  //! TransposedMatrix-Vector multiplication for arbitrary vectors that support operator[]
  template <class XX, class YY>
  inline std::enable_if_t<XT::Common::VectorAbstraction<XX>::is_vector && XT::Common::VectorAbstraction<YY>::is_vector,
                          void>
  mtv(const XX& xx, YY& yy) const
  {
    std::fill(yy.begin(), yy.end(), ScalarType(0));
    const auto& entries = *entries_;
    const auto& row_pointers = *row_pointers_;
    const auto& column_indices = *column_indices_;
    for (size_t rr = 0; rr < num_rows_; ++rr) {
      const size_t end = row_pointers[rr + 1];
      for (size_t kk = row_pointers[rr]; kk < end; ++kk)
        yy[column_indices[kk]] += entries[kk] * xx[rr];
    }
  }

  inline void add_to_entry(const size_t rr, const size_t cc, const ScalarType& value)
  {
    ensure_uniqueness();
    internal::LockGuard DUNE_UNUSED(lock)(mutexes_, rr);
    entries_->operator[](get_entry_index(rr, cc)) += value;
  }

  inline ScalarType get_entry(const size_t rr, const size_t cc) const
  {
    const size_t index = get_entry_index(rr, cc, false);
    return index == size_t(-1) ? ScalarType(0) : entries_->operator[](index);
  }

  inline void set_entry(const size_t rr, const size_t cc, const ScalarType value)
  {
    ensure_uniqueness();
    entries_->operator[](get_entry_index(rr, cc)) = value;
  }

  inline void clear_row(const size_t rr)
  {
    ensure_uniqueness();
    std::fill(entries_->begin() + row_pointers_->operator[](rr),
              entries_->begin() + row_pointers_->operator[](rr + 1),
              ScalarType(0));
  }

  inline void clear_col(const size_t cc)
  {
    ensure_uniqueness();
    for (size_t kk = 0; kk < entries_->size(); ++kk) {
      if (column_indices_->operator[](kk) == cc)
        entries_->operator[](kk) = ScalarType(0);
    }
  }

  inline void unit_row(const size_t rr)
  {
    clear_row(rr);
    set_entry(rr, rr, ScalarType(1));
  }

  inline void unit_col(const size_t cc)
  {
    clear_col(cc);
    set_entry(cc, cc, ScalarType(1));
  }

  bool valid() const
  {
    // iterate over non-zero entries
    for (const auto& entry : *entries_)
      if (XT::Common::isnan(std::real(entry)) || XT::Common::isnan(std::imag(entry))
          || XT::Common::isinf(std::abs(entry)))
        return false;
    return true;
  }

  virtual size_t non_zeros() const override final
  {
    return entries_->size();
  }

  virtual SparsityPatternDefault pattern(const bool prune = false,
                                         const EpsType eps = Common::FloatCmp::DefaultEpsilon<ScalarType>::value()
                                                             / 1000.) const override
  {
    SparsityPatternDefault ret(num_rows_);
    for (size_t rr = 0; rr < num_rows_; ++rr) {
      for (size_t kk = row_pointers_->operator[](rr); kk < row_pointers_->operator[](rr + 1); ++kk) {
        if (!prune || Common::FloatCmp::ne(entries_->operator[](kk), ScalarType(0), 0., eps / num_cols_))
          ret.insert(rr, column_indices_->operator[](kk));
      }
    }
    return ret;
  } // ... pattern(...)

  /// \}

  template <class MatrixType>
  void rightmultiply(const MatrixType& other)
  {
    ensure_uniqueness();
    auto new_entries = std::make_shared<EntriesVectorType>();
    auto new_row_pointers = std::make_shared<IndexVectorType>(num_rows_ + 1);
    auto new_column_indices = std::make_shared<IndexVectorType>();
    new_entries->reserve(entries_->size());
    new_column_indices->reserve(column_indices_->size());
    ScalarType new_entry(0);
    size_t index = 0;
    for (size_t rr = 0; rr < num_rows_; ++rr) {
      (*new_row_pointers)[rr] = index;
      for (size_t cc = 0; cc < num_cols_; ++cc) {
        new_entry = 0;
        for (size_t kk = (*row_pointers_)[rr]; kk < (*row_pointers_)[rr + 1]; ++kk)
          new_entry += (*entries_)[kk] * other[(*column_indices_)[kk]][cc];
        if (XT::Common::FloatCmp::ne(new_entry, 0., 0., eps_ / num_cols_)) {
          new_entries->push_back(new_entry);
          new_column_indices->push_back(cc);
          ++index;
        }
      } // cc
    } // rr
    (*new_row_pointers)[num_rows_] = new_column_indices->size();
    entries_ = new_entries;
    row_pointers_ = new_row_pointers;
    column_indices_ = new_column_indices;
  } // void rightmultiply(...)

  using MatrixInterfaceType::operator+;
  using MatrixInterfaceType::operator-;
  using MatrixInterfaceType::operator+=;
  using MatrixInterfaceType::operator-=;

  const EntriesVectorType& entries() const
  {
    return *entries_;
  }

  const IndexVectorType& row_pointers() const
  {
    return *row_pointers_;
  }

  const IndexVectorType& column_indices() const
  {
    return *column_indices_;
  }

  const EpsType& eps()
  {
    return eps_;
  }

protected:
  inline void ensure_uniqueness() const
  {
    if (!entries_.unique()) {
      assert(!unshareable_);
      const internal::VectorLockGuard DUNE_UNUSED(guard)(mutexes_);
      if (!entries_.unique()) {
        entries_ = std::make_shared<EntriesVectorType>(*entries_);
        row_pointers_ = std::make_shared<IndexVectorType>(*row_pointers_);
        column_indices_ = std::make_shared<IndexVectorType>(*column_indices_);
        mutexes_ = mutexes_ ? std::make_shared<std::vector<std::mutex>>(mutexes_->size()) : nullptr;
      }
    }
  } // ... ensure_uniqueness(...)

private:
  size_t get_entry_index(const size_t rr, const size_t cc, const bool throw_if_not_in_pattern = true) const
  {
    const auto& row_offset = row_pointers_->operator[](rr);
    const auto& next_row_offset = row_pointers_->operator[](rr + 1);
    const auto column_indices_it = column_indices_->begin() + row_offset;
    const auto column_indices_it_end = column_indices_->begin() + next_row_offset;
    const auto entry_it = std::lower_bound(column_indices_it, column_indices_it_end, cc);
    if (entry_it != column_indices_it_end && *entry_it == cc)
      return row_offset + std::distance(column_indices_it, entry_it);
    if (throw_if_not_in_pattern)
      DUNE_THROW(Common::Exceptions::index_out_of_range, "Entry is not in the sparsity pattern!");
    return size_t(-1);
  }

  size_t num_rows_, num_cols_;
  mutable std::shared_ptr<EntriesVectorType> entries_;
  mutable std::shared_ptr<IndexVectorType> row_pointers_;
  mutable std::shared_ptr<IndexVectorType> column_indices_;
  mutable std::shared_ptr<std::vector<std::mutex>> mutexes_;
  mutable EpsType eps_;
  mutable bool unshareable_;
}; // class CommonSparseMatrix

/**
 * \brief A sparse matrix implementation of the MatrixInterface with column major memory layout.
 */
template <class ScalarImp>
class CommonSparseMatrix<ScalarImp, SparseFormat::csc>
    : public MatrixInterface<internal::CommonSparseMatrixTraits<ScalarImp, SparseFormat::csc>, ScalarImp>
{
  typedef CommonSparseMatrix<ScalarImp, SparseFormat::csc> ThisType;
  typedef MatrixInterface<internal::CommonSparseMatrixTraits<ScalarImp, SparseFormat::csc>, ScalarImp>
      MatrixInterfaceType;

public:
  typedef internal::CommonSparseMatrixTraits<ScalarImp, SparseFormat::csc> Traits;
  typedef typename Traits::EntriesVectorType EntriesVectorType;
  typedef typename Traits::ScalarType ScalarType;
  typedef typename Traits::RealType RealType;
  typedef typename Traits::EpsType EpsType;
  typedef std::vector<size_t> IndexVectorType;

  /**
  * \brief This is the constructor of interest which creates a sparse matrix.
  */
  CommonSparseMatrix(const size_t rr,
                     const size_t cc,
                     const SparsityPatternDefault& patt,
                     const size_t num_mutexes = 1,
                     const EpsType eps = Common::FloatCmp::DefaultEpsilon<ScalarType>::value() / 1000.)
    : num_rows_(rr)
    , num_cols_(cc)
    , entries_(std::make_shared<EntriesVectorType>())
    , column_pointers_(std::make_shared<IndexVectorType>(num_cols_ + 1, 0))
    , row_indices_(std::make_shared<IndexVectorType>())
    , mutexes_(num_mutexes > 0 ? std::make_shared<std::vector<std::mutex>>(num_mutexes) : nullptr)
    , eps_(eps)
    , unshareable_(false)
  {
    if (num_rows_ > 0 && num_cols_ > 0) {
      if (size_t(patt.size()) != num_rows_)
        DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                   "The size of the pattern (" << patt.size() << ") does not match the number of rows of this ("
                                               << num_rows_
                                               << ")!");
      for (size_t col = 0; col < num_cols_; ++col) {
        for (size_t row = 0; row < num_rows_; ++row) {
          const auto& column_indices = patt.inner(row);
          if (std::find(column_indices.begin(), column_indices.end(), col) != column_indices.end())
            row_indices_->push_back(row);
        } // row
        (*column_pointers_)[col + 1] = row_indices_->size();
      } // col
      entries_->resize(row_indices_->size());
    }
  } // CommonSparseMatrix(rr, cc, patt, num_mutexes)

  CommonSparseMatrix(const size_t rr = 0,
                     const size_t cc = 0,
                     const ScalarType& value = ScalarType(0),
                     const size_t num_mutexes = 1,
                     const EpsType eps = Common::FloatCmp::DefaultEpsilon<ScalarType>::value() / 1000.)
    : num_rows_(rr)
    , num_cols_(cc)
    , entries_(std::make_shared<EntriesVectorType>(
          XT::Common::FloatCmp::eq(value, ScalarType(0.), 0., eps / num_cols_) ? 0 : num_rows_ * num_cols_, value))
    , column_pointers_(std::make_shared<IndexVectorType>(num_cols_ + 1, 0))
    , row_indices_(std::make_shared<IndexVectorType>())
    , mutexes_(num_mutexes > 0 ? std::make_shared<std::vector<std::mutex>>(num_mutexes) : nullptr)
    , eps_(eps)
    , unshareable_(false)
  {
    if (XT::Common::FloatCmp::ne(value, ScalarType(0.), 0., eps / num_cols_)) {
      IndexVectorType column_row_indices(num_rows_);
      for (size_t row = 0; row < num_rows_; ++row)
        column_row_indices[row] = row;
      for (size_t col = 0; col < num_cols_; ++col) {
        column_pointers_->operator[](col + 1) = (col + 1) * num_rows_;
        row_indices_->insert(row_indices_->end(), column_row_indices.begin(), column_row_indices.end());
      }
    }
  }

  CommonSparseMatrix(const size_t rr, const size_t cc, const size_t num_mutexes)
    : num_rows_(rr)
    , num_cols_(cc)
    , entries_(std::make_shared<EntriesVectorType>())
    , column_pointers_(std::make_shared<IndexVectorType>(num_cols_ + 1))
    , row_indices_(std::make_shared<IndexVectorType>())
    , mutexes_(num_mutexes > 0 ? std::make_shared<std::vector<std::mutex>>(num_mutexes) : nullptr)
    , unshareable_(false)
  {
  }

  CommonSparseMatrix(const ThisType& other)
    : num_rows_(other.num_rows_)
    , num_cols_(other.num_cols_)
    , entries_(other.unshareable_ ? std::make_shared<EntriesVectorType>(*other.entries_) : other.entries_)
    , column_pointers_(other.unshareable_ ? std::make_shared<IndexVectorType>(*other.column_pointers_)
                                          : other.column_pointers_)
    , row_indices_(other.unshareable_ ? std::make_shared<IndexVectorType>(*other.row_indices_) : other.row_indices_)
    , mutexes_(other.unshareable_ ? std::make_shared<std::vector<std::mutex>>(other.mutexes_->size()) : other.mutexes_)
    , unshareable_(false)
  {
  }

  template <class OtherMatrixType>
  explicit CommonSparseMatrix(
      const OtherMatrixType& mat,
      const typename std::enable_if<Common::MatrixAbstraction<OtherMatrixType>::is_matrix, bool>::type prune = false,
      const EpsType eps = Common::FloatCmp::DefaultEpsilon<ScalarType>::value() / 1000.,
      const size_t num_mutexes = 1)
    : num_rows_(Common::MatrixAbstraction<OtherMatrixType>::rows(mat))
    , num_cols_(Common::MatrixAbstraction<OtherMatrixType>::cols(mat))
    , entries_(std::make_shared<EntriesVectorType>())
    , column_pointers_(std::make_shared<IndexVectorType>(num_cols_ + 1, 0))
    , row_indices_(std::make_shared<IndexVectorType>())
    , mutexes_(num_mutexes > 0 ? std::make_shared<std::vector<std::mutex>>(num_mutexes) : nullptr)
    , eps_(eps)
    , unshareable_(false)
  {
    entries_->reserve(0.1 * num_rows_ * num_cols_);
    row_indices_->reserve(0.1 * num_rows_ * num_cols_);
    size_t index = 0;
    for (size_t cc = 0; cc < num_cols_; ++cc) {
      for (size_t rr = 0; rr < num_rows_; ++rr) {
        const auto& value = Common::MatrixAbstraction<OtherMatrixType>::get_entry(mat, rr, cc);
        if (!prune || XT::Common::FloatCmp::ne(
                          value, ScalarType(0), 0., eps / Common::MatrixAbstraction<OtherMatrixType>::cols(mat))) {
          entries_->push_back(value);
          row_indices_->push_back(rr);
          ++index;
        }
      } // rr
      (*column_pointers_)[cc + 1] = index;
    } // cc
  } // CommonSparseMatrix(...)

  template <class OtherMatrixType>
  explicit CommonSparseMatrix(
      const OtherMatrixType& mat,
      const typename std::enable_if<Common::MatrixAbstraction<OtherMatrixType>::is_matrix, bool>::type prune,
      const size_t num_mutexes)
    : CommonSparseMatrix(mat,
                         prune,
                         Common::FloatCmp::DefaultEpsilon<ScalarType>::value()
                             / (1000. * Common::MatrixAbstraction<OtherMatrixType>::cols(mat)),
                         num_mutexes)
  {
  } // CommonSparseMatrix(...)

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      num_rows_ = other.num_rows_;
      num_cols_ = other.num_cols_;
      entries_ = other.unshareable_ ? std::make_shared<EntriesVectorType>(*other.entries_) : other.entries_;
      column_pointers_ = other.column_pointers_;
      row_indices_ = other.row_indices_;
      mutexes_ = other.mutexes_
                     ? (other.unshareable_ ? std::make_shared<std::vector<std::mutex>>(other.mutexes_->size())
                                           : other.mutexes_)
                     : nullptr;
      eps_(other.eps_);
      unshareable_ = false;
    }
    return *this;
  }

  template <class OtherMatrixImp>
  ThisType& operator=(const OtherMatrixImp& other)
  {
    typedef Common::MatrixAbstraction<OtherMatrixImp> MatAbstrType;
    clear();
    num_rows_ = MatAbstrType::rows(other);
    num_cols_ = MatAbstrType::cols(other);
    for (size_t cc = 0; cc < num_cols_; ++cc) {
      for (size_t rr = 0; rr < num_rows_; ++rr) {
        if (XT::Common::FloatCmp::ne(MatAbstrType::get_entry(other, rr, cc), ScalarType(0), 0., eps_ / num_cols_)) {
          entries_->push_back(MatAbstrType::get_entry(other, rr, cc));
          row_indices_->push_back(rr);
        }
      } // rr
      (*column_pointers_)[cc + 1] = row_indices_->size();
    } // cc
    return *this;
  }

  void deep_copy(const ThisType& other)
  {
    ensure_uniqueness();
    num_rows_ = other.num_rows_;
    num_cols_ = other.num_cols_;
    *entries_ = *other.entries_;
    *column_pointers_ = *other.column_pointers_;
    *row_indices_ = *other.row_indices_;
  }

  void clear()
  {
    ensure_uniqueness();
    entries_->clear();
    std::fill(column_pointers_->begin(), column_pointers_->end(), 0);
    row_indices_->clear();
  }

  /// \name Required by ContainerInterface.
  /// \{
  inline ThisType copy() const
  {
    ThisType ret(*this);
    ret.entries_ = std::make_shared<EntriesVectorType>(*entries_);
    ret.column_pointers_ = std::make_shared<IndexVectorType>(*column_pointers_);
    ret.row_indices_ = std::make_shared<IndexVectorType>(*row_indices_);
    ret.mutexes_ = mutexes_ ? std::make_shared<std::vector<std::mutex>>(mutexes_->size()) : nullptr;
    return ret;
  }

  inline void scal(const ScalarType& alpha)
  {
    ensure_uniqueness();
    const internal::VectorLockGuard DUNE_UNUSED(guard)(mutexes_);
    std::transform(
        entries_->begin(), entries_->end(), entries_->begin(), std::bind1st(std::multiplies<ScalarType>(), alpha));
  }

  inline void axpy(const ScalarType& alpha, const ThisType& xx)
  {
    ensure_uniqueness();
    const internal::VectorLockGuard DUNE_UNUSED(guard)(mutexes_);
    assert(has_equal_shape(xx));
    const auto& xx_entries = *xx.entries_;
    for (size_t ii = 0; ii < entries_->size(); ++ii)
      entries_->operator[](ii) += alpha * xx_entries[ii];
  }

  inline bool has_equal_shape(const ThisType& other) const
  {
    return (rows() == other.rows()) && (cols() == other.cols());
  }

  /// \}
  /// \name Required by MatrixInterface.
  /// \{

  inline size_t rows() const
  {
    return num_rows_;
  }

  inline size_t cols() const
  {
    return num_cols_;
  }

  //! Matrix-Vector multiplication for arbitrary vectors that support operator[]
  template <class XX, class YY>
  inline std::enable_if_t<!std::is_base_of<CommonSparseVector<ScalarType>, XX>::value
                              && !std::is_base_of<CommonSparseVector<ScalarType>, YY>::value
                              && XT::Common::VectorAbstraction<XX>::is_vector
                              && XT::Common::VectorAbstraction<YY>::is_vector,
                          void>
  mv(const XX& xx, YY& yy) const
  {
    std::fill(yy.begin(), yy.end(), ScalarType(0));
    const auto& entries = *entries_;
    const auto& column_pointers = *column_pointers_;
    const auto& row_indices = *row_indices_;
    for (size_t cc = 0; cc < num_cols_; ++cc) {
      const size_t end = column_pointers[cc + 1];
      for (size_t kk = column_pointers[cc]; kk < end; ++kk)
        yy[row_indices[kk]] += entries[kk] * xx[cc];
    }
  }

  void mv(const CommonSparseVector<ScalarType>& xx, CommonSparseVector<ScalarType>& yy) const
  {
    yy.clear();
    const auto& entries = *entries_;
    const auto& column_pointers = *column_pointers_;
    const auto& row_indices = *row_indices_;
    const auto& vec_entries = xx.entries();
    const auto& vec_indices = xx.indices();
    thread_local std::vector<ScalarType> tmp_vec;
    tmp_vec.resize(num_rows_);
    std::fill(tmp_vec.begin(), tmp_vec.end(), 0.);
    for (size_t ii = 0; ii < vec_entries.size(); ++ii) {
      const size_t cc = vec_indices[ii];
      const size_t end = column_pointers[cc + 1];
      for (size_t kk = column_pointers[cc]; kk < end; ++kk)
        tmp_vec[row_indices[kk]] += entries[kk] * vec_entries[ii];
    }
    for (size_t cc = 0; cc < num_cols_; ++cc)
      if (XT::Common::FloatCmp::ne(tmp_vec[cc], ScalarType(0.), 0., eps_))
        yy.set_new_entry(cc, tmp_vec[cc]);
  }

  //! TransposedMatrix-Vector multiplication for arbitrary vectors that support operator[]
  template <class XX, class YY>
  inline std::enable_if_t<!std::is_base_of<CommonSparseVector<ScalarType>, XX>::value
                              && !std::is_base_of<CommonSparseVector<ScalarType>, YY>::value
                              && XT::Common::VectorAbstraction<XX>::is_vector
                              && XT::Common::VectorAbstraction<YY>::is_vector,
                          void>
  mtv(const XX& xx, YY& yy) const
  {
    std::fill(yy.begin(), yy.end(), ScalarType(0));
    const auto& entries = *entries_;
    const auto& column_pointers = *column_pointers_;
    const auto& row_indices = *row_indices_;
    for (size_t cc = 0; cc < num_cols_; ++cc) {
      const size_t end = column_pointers[cc + 1];
      for (size_t kk = column_pointers[cc]; kk < end; ++kk)
        yy[cc] += entries[kk] * xx[row_indices[kk]];
    }
  }

  void mtv(const CommonSparseVector<ScalarType>& xx, CommonSparseVector<ScalarType>& yy) const
  {
    yy.clear();
    const auto& entries = *entries_;
    const auto& column_pointers = *column_pointers_;
    const auto& row_indices = *row_indices_;
    for (size_t cc = 0; cc < num_cols_; ++cc) {
      const size_t end = column_pointers[cc + 1];
      ScalarType new_entry(0);
      for (size_t kk = column_pointers[cc]; kk < end; ++kk)
        new_entry += entries[kk] * xx.get_entry(row_indices[kk]);
      if (XT::Common::FloatCmp::ne(new_entry, ScalarType(0.), 0., eps_))
        yy.set_new_entry(cc, new_entry);
    }
  }

  inline void add_to_entry(const size_t rr, const size_t cc, const ScalarType& value)
  {
    ensure_uniqueness();
    internal::LockGuard DUNE_UNUSED(lock)(mutexes_, rr);
    entries_->operator[](get_entry_index(rr, cc)) += value;
  }

  inline ScalarType get_entry(const size_t rr, const size_t cc) const
  {
    const size_t index = get_entry_index(rr, cc, false);
    return index == size_t(-1) ? ScalarType(0) : (*entries_)[index];
  }

  inline void set_entry(const size_t rr, const size_t cc, const ScalarType value)
  {
    ensure_uniqueness();
    entries_->operator[](get_entry_index(rr, cc)) = value;
  }

  inline void clear_row(const size_t rr)
  {
    ensure_uniqueness();
    for (size_t kk = 0; kk < entries_->size(); ++kk) {
      if ((*row_indices_)[kk] == rr)
        (*entries_)[kk] = ScalarType(0);
    }
  }

  inline void clear_col(const size_t cc)
  {
    ensure_uniqueness();
    std::fill(
        entries_->begin() + (*column_pointers_)[cc], entries_->begin() + (*column_pointers_)[cc + 1], ScalarType(0));
  }

  inline void unit_row(const size_t rr)
  {
    clear_row(rr);
    set_entry(rr, rr, ScalarType(1));
  }

  inline void unit_col(const size_t cc)
  {
    clear_col(cc);
    set_entry(cc, cc, ScalarType(1));
  }

  bool valid() const
  {
    // iterate over non-zero entries
    for (const auto& entry : *entries_)
      if (XT::Common::isnan(std::real(entry)) || XT::Common::isnan(std::imag(entry))
          || XT::Common::isinf(std::abs(entry)))
        return false;
    return true;
  }

  virtual size_t non_zeros() const override final
  {
    return entries_->size();
  }

  virtual SparsityPatternDefault pattern(const bool prune = false,
                                         const EpsType eps = Common::FloatCmp::DefaultEpsilon<ScalarType>::value()
                                                             / 1000.) const override
  {
    SparsityPatternDefault ret(num_rows_);
    for (size_t cc = 0; cc < num_cols_; ++cc) {
      for (size_t kk = (*column_pointers_)[cc]; kk < (*column_pointers_)[cc + 1]; ++kk) {
        if (!prune || Common::FloatCmp::ne((*entries_)[kk], ScalarType(0), 0., eps / num_cols_))
          ret.insert((*row_indices_)[kk], cc);
      }
    } // cc
    return ret;
  } // ... pattern(...)

  /// \}

  template <class OtherMatrixImp>
  typename std::enable_if_t<XT::Common::MatrixAbstraction<OtherMatrixImp>::is_matrix
                                && !(std::is_base_of<ThisType, OtherMatrixImp>::value),
                            void>
  rightmultiply(const OtherMatrixImp& other)
  {
    ensure_uniqueness();
    EntriesVectorType new_entries;
    new_entries.reserve(entries_->size());
    IndexVectorType new_column_pointers(num_cols_ + 1, 0);
    IndexVectorType new_row_indices;
    new_row_indices.reserve(row_indices_->size());
    ScalarType new_entry(0);
    size_t index = 0;
    for (size_t cc = 0; cc < num_cols_; ++cc) {
      for (size_t rr = 0; rr < num_rows_; ++rr) {
        new_entry = 0;
        for (size_t col = 0; col < num_cols_; ++col) {
          for (size_t kk = (*column_pointers_)[col]; kk < (*column_pointers_)[col + 1]; ++kk)
            if ((*row_indices_)[kk] == rr)
              new_entry += (*entries_)[kk] * XT::Common::MatrixAbstraction<OtherMatrixImp>::get_entry(other, col, cc);
        } // col
        if (XT::Common::FloatCmp::ne(new_entry, 0., 0., eps_ / num_cols_)) {
          new_entries.push_back(new_entry);
          new_row_indices.push_back(rr);
          ++index;
        }
      } // rr
      new_column_pointers[cc + 1] = index;
    } // cc
    *entries_ = new_entries;
    *column_pointers_ = new_column_pointers;
    *row_indices_ = new_row_indices;
  } // void rightmultiply(...)

  void rightmultiply(const ThisType& other)
  {
    ensure_uniqueness();
    thread_local EntriesVectorType new_entries;
    thread_local IndexVectorType new_column_pointers(num_cols_ + 1, 0);
    thread_local IndexVectorType new_row_indices;
    new_entries.clear();
    new_row_indices.clear();
    const auto& other_entries = other.entries();
    const auto& other_column_pointers = other.column_pointers();
    const auto& other_row_indices = other.row_indices();
    thread_local std::vector<ScalarType> dense_column(num_rows_, 0.);
    dense_column.resize(num_rows_);
    for (size_t cc = 0; cc < other.cols(); ++cc) {
      std::fill(dense_column.begin(), dense_column.end(), 0.);
      size_t other_col_end = other_column_pointers[cc + 1];
      for (size_t kk = other_column_pointers[cc]; kk < other_col_end; ++kk) {
        size_t col = other_row_indices[kk];
        size_t col_end = (*column_pointers_)[col + 1];
        for (size_t ll = (*column_pointers_)[col]; ll < col_end; ++ll)
          dense_column[(*row_indices_)[ll]] += (*entries_)[ll] * other_entries[kk];
      } // kk
      for (size_t rr = 0; rr < num_rows_; ++rr) {
        if (XT::Common::FloatCmp::ne(dense_column[rr], ScalarType(0.), 0., eps_ / num_cols_)) {
          new_entries.push_back(dense_column[rr]);
          new_row_indices.push_back(rr);
        }
      } // ii
      new_column_pointers[cc + 1] = new_row_indices.size();
    } // cc
    *entries_ = new_entries;
    *column_pointers_ = new_column_pointers;
    *row_indices_ = new_row_indices;
  } // void rightmultiply(...)

  using MatrixInterfaceType::operator+;
  using MatrixInterfaceType::operator-;
  using MatrixInterfaceType::operator+=;
  using MatrixInterfaceType::operator-=;

  EntriesVectorType& entries()
  {
    return *entries_;
  }

  const EntriesVectorType& entries() const
  {
    return *entries_;
  }

  const IndexVectorType& column_pointers() const
  {
    return *column_pointers_;
  }

  IndexVectorType& column_pointers()
  {
    return *column_pointers_;
  }

  IndexVectorType& row_indices()
  {
    return *row_indices_;
  }

  const IndexVectorType& row_indices() const
  {
    return *row_indices_;
  }

  const EpsType& eps()
  {
    return eps_;
  }

protected:
  inline void ensure_uniqueness() const
  {
    if (!entries_.unique()) {
      assert(!unshareable_);
      const internal::VectorLockGuard DUNE_UNUSED(guard)(mutexes_);
      if (!entries_.unique()) {
        entries_ = std::make_shared<EntriesVectorType>(*entries_);
        column_pointers_ = std::make_shared<IndexVectorType>(*column_pointers_);
        row_indices_ = std::make_shared<IndexVectorType>(*row_indices_);
        mutexes_ = mutexes_ ? std::make_shared<std::vector<std::mutex>>(mutexes_->size()) : nullptr;
      }
    }
  } // ... ensure_uniqueness(...)

private:
  size_t get_entry_index(const size_t rr, const size_t cc, const bool throw_if_not_in_pattern = true) const
  {
    const auto& column_offset = column_pointers_->operator[](cc);
    const auto& next_column_offset = column_pointers_->operator[](cc + 1);
    const auto it = row_indices_->begin() + column_offset;
    const auto it_end = row_indices_->begin() + next_column_offset;
    const auto entry_it = std::lower_bound(it, it_end, rr);
    if (entry_it != it_end && *entry_it == rr)
      return column_offset + std::distance(it, entry_it);
    if (throw_if_not_in_pattern)
      DUNE_THROW(Common::Exceptions::index_out_of_range, "Entry is not in the sparsity pattern!");
    return size_t(-1);
  }

  size_t num_rows_, num_cols_;
  mutable std::shared_ptr<EntriesVectorType> entries_;
  mutable std::shared_ptr<IndexVectorType> column_pointers_;
  mutable std::shared_ptr<IndexVectorType> row_indices_;
  mutable std::shared_ptr<std::vector<std::mutex>> mutexes_;
  mutable EpsType eps_;
  mutable bool unshareable_;
}; // class CommonSparseMatrix<..., SparseFormat::csc>

/**
 * \brief A matrix implementation checking whether the matrix is sparse enough to use sparse matrix operations.
 */
template <class DenseMatrixImp, class SparseMatrixImp>
class CommonSparseOrDenseMatrix
    : public MatrixInterface<internal::CommonSparseOrDenseMatrixTraits<DenseMatrixImp, SparseMatrixImp>,
                             typename SparseMatrixImp::ScalarType>
{
  typedef CommonSparseOrDenseMatrix<DenseMatrixImp, SparseMatrixImp> ThisType;
  typedef MatrixInterface<internal::CommonSparseOrDenseMatrixTraits<DenseMatrixImp, SparseMatrixImp>,
                          typename SparseMatrixImp::ScalarType>
      MatrixInterfaceType;

public:
  typedef internal::CommonSparseOrDenseMatrixTraits<DenseMatrixImp, SparseMatrixImp> Traits;
  typedef typename Traits::DenseMatrixType DenseMatrixType;
  typedef typename Traits::SparseMatrixType SparseMatrixType;
  typedef typename Traits::ScalarType ScalarType;
  typedef typename Traits::RealType RealType;
  typedef typename Traits::EpsType EpsType;
  static constexpr double sparse_limit = 0.1;

  /**
  * \brief This is the constructor of interest which creates a sparse matrix.
  */
  CommonSparseOrDenseMatrix(const size_t rr,
                            const size_t cc,
                            const SparsityPatternDefault& patt,
                            const size_t num_mutexes = 1,
                            EpsType eps = Common::FloatCmp::DefaultEpsilon<ScalarType>::value() / 1000.)
    : num_rows_(rr)
    , num_cols_(cc)
  {
    size_t nnz = 0;
    for (size_t row = 0; row < num_rows_; ++row)
      nnz += patt.inner(row).size();
    size_t num_entries = rr * cc;
    double density = double(nnz) / double(num_entries);
    sparse_ = density < sparse_limit;
    if (sparse_) {
      sparse_matrix_ = SparseMatrixType(rr, cc, patt, num_mutexes, eps);
      dense_matrix_ = DenseMatrixType(0, 0, patt, num_mutexes);
    } else {
      sparse_matrix_ = SparseMatrixType(0, 0, patt, num_mutexes, eps);
      dense_matrix_ = DenseMatrixType(rr, cc, patt, num_mutexes);
    }
  } // CommonSparseOrDenseMatrix(rr, cc, patt, num_mutexes)

  CommonSparseOrDenseMatrix(const size_t rr = 0,
                            const size_t cc = 0,
                            const ScalarType& value = ScalarType(0),
                            const size_t num_mutexes = 1,
                            bool use_sparse_if_zero = true,
                            const EpsType eps = Common::FloatCmp::DefaultEpsilon<ScalarType>::value() / 1000.)
    : num_rows_(rr)
    , num_cols_(cc)
  {
    if (XT::Common::FloatCmp::ne(value, ScalarType(0.), 0., eps / num_cols_) || !use_sparse_if_zero) {
      sparse_matrix_ = SparseMatrixType(0, 0, value, num_mutexes, eps);
      dense_matrix_ = DenseMatrixType(rr, cc, value, num_mutexes);
    } else {
      sparse_matrix_ = SparseMatrixType(rr, cc, value, num_mutexes, eps);
      dense_matrix_ = DenseMatrixType(0, 0, value, num_mutexes);
      sparse_ = true;
    }
  }

  CommonSparseOrDenseMatrix(const size_t rr, const size_t cc, const size_t num_mutexes, bool use_sparse = true)
    : CommonSparseOrDenseMatrix(rr, cc, ScalarType(0.), num_mutexes, use_sparse)
  {
  }

  CommonSparseOrDenseMatrix(const ThisType& other)
    : num_rows_(other.num_rows_)
    , num_cols_(other.num_cols_)
    , sparse_(other.sparse_)
    , sparse_matrix_(other.sparse_matrix_)
    , dense_matrix_(other.dense_matrix_)
  {
  }

  template <class OtherMatrixType>
  explicit CommonSparseOrDenseMatrix(
      const OtherMatrixType& mat,
      const typename std::enable_if<Common::MatrixAbstraction<OtherMatrixType>::is_matrix, bool>::type prune = false,
      const EpsType eps_in = Common::FloatCmp::DefaultEpsilon<ScalarType>::value() / 1000.,
      const size_t num_mutexes = 1,
      const bool prune_dense = false)
    : num_rows_(Common::MatrixAbstraction<OtherMatrixType>::rows(mat))
    , num_cols_(Common::MatrixAbstraction<OtherMatrixType>::cols(mat))
  {
    auto eps = eps_in / Common::MatrixAbstraction<OtherMatrixType>::cols(mat);
    // check sparsity (if prune = false, the sparsity may not be checked correctly)
    size_t nnz = 0.;
    for (size_t rr = 0; rr < num_rows_; ++rr)
      for (size_t cc = 0; cc < num_cols_; ++cc) {
        const auto& value = Common::MatrixAbstraction<OtherMatrixType>::get_entry(mat, rr, cc);
        nnz += XT::Common::FloatCmp::ne(value, ScalarType(0), 0., eps);
      }
    double density = double(nnz) / double(num_rows_ * num_cols_);
    sparse_ = density < sparse_limit;
    if (sparse_) {
      sparse_matrix_ = SparseMatrixType(mat, prune, eps_in, num_mutexes);
      dense_matrix_ = DenseMatrixType(0, 0, num_mutexes);
    } else {
      sparse_matrix_ = SparseMatrixType(0, 0, num_mutexes);
      dense_matrix_ = DenseMatrixType(mat, prune_dense ? prune : false, eps_in, num_mutexes);
    } // else (sparse_)
  } // CommonSparseOrDenseMatrix(...)

  template <class OtherMatrixType>
  explicit CommonSparseOrDenseMatrix(
      const OtherMatrixType& mat,
      const typename std::enable_if<Common::MatrixAbstraction<OtherMatrixType>::is_matrix, bool>::type prune,
      const size_t num_mutexes)
    : CommonSparseOrDenseMatrix(mat,
                                prune,
                                Common::FloatCmp::DefaultEpsilon<ScalarType>::value
                                    / (10. * Common::MatrixAbstraction<OtherMatrixType>::cols(mat)),
                                num_mutexes)
  {
  } // CommonSparseOrDenseMatrix(...)

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      num_rows_ = other.num_rows_;
      num_cols_ = other.num_cols_;
      sparse_ = other.sparse_;
      sparse_matrix_ = other.sparse_matrix_;
      dense_matrix_ = other.dense_matrix_;
    }
    return *this;
  }

  ThisType& operator=(const DenseMatrixType& other)
  {
    size_t nnz = 0.;
    auto tol = sparse_matrix_.eps() / Common::MatrixAbstraction<DenseMatrixType>::cols(other);
    for (size_t rr = 0; rr < num_rows_; ++rr)
      for (size_t cc = 0; cc < num_cols_; ++cc) {
        const auto& value = Common::MatrixAbstraction<DenseMatrixType>::get_entry(other, rr, cc);
        nnz += XT::Common::FloatCmp::ne(value, ScalarType(0), 0., tol);
      }
    double density = double(nnz) / (Common::MatrixAbstraction<DenseMatrixType>::rows(other)
                                    * Common::MatrixAbstraction<DenseMatrixType>::cols(other));
    sparse_ = density < sparse_limit;
    if (sparse_)
      sparse_matrix_ = other;
    else
      dense_matrix_ = other;
    return *this;
  }

  void deep_copy(const ThisType& other)
  {
    num_rows_ = other.num_rows_;
    num_cols_ = other.num_cols_;
    sparse_ = other.sparse_;
    sparse_matrix_.deep_copy(other.sparse_matrix_);
    dense_matrix_.deep_copy(other.dense_matrix_);
  }

  //  void clear()
  //  {
  //    ensure_uniqueness();
  //    entries_->clear();
  //    std::fill(column_pointers_->begin(), column_pointers_->end(), 0);
  //    row_indices_->clear();
  //  }

  /// \name Required by ContainerInterface.
  /// \{
  inline ThisType copy() const
  {
    ThisType ret(*this);
    ret.sparse_matrix_ = sparse_matrix_.copy();
    ret.dense_matrix_ = dense_matrix_.copy();
    return ret;
  }

  inline void scal(const ScalarType& alpha)
  {
    sparse_ ? sparse_matrix_.scal(alpha) : dense_matrix_.scal(alpha);
  }

  inline void axpy(const ScalarType& alpha, const ThisType& xx)
  {
    assert(sparse_ == xx.sparse_);
    sparse_ ? sparse_matrix_.axpy(alpha, xx.sparse_matrix_) : dense_matrix_.axpy(alpha, xx.dense_matrix_);
  }

  inline bool has_equal_shape(const ThisType& other) const
  {
    return (rows() == other.rows()) && (cols() == other.cols());
  }

  /// \}
  /// \name Required by MatrixInterface.
  /// \{

  inline size_t rows() const
  {
    return num_rows_;
  }

  inline size_t cols() const
  {
    return num_cols_;
  }

  //! Matrix-Vector multiplication for arbitrary vectors that support operator[]
  template <class XX, class YY>
  inline std::enable_if_t<XT::Common::VectorAbstraction<XX>::is_vector && XT::Common::VectorAbstraction<YY>::is_vector,
                          void>
  mv(const XX& xx, YY& yy) const
  {
    sparse_ ? sparse_matrix_.mv(xx, yy) : dense_matrix_.mv(xx, yy);
  }

  //! TransposedMatrix-Vector multiplication for arbitrary vectors that support operator[]
  template <class XX, class YY>
  inline std::enable_if_t<XT::Common::VectorAbstraction<XX>::is_vector && XT::Common::VectorAbstraction<YY>::is_vector,
                          void>
  mtv(const XX& xx, YY& yy) const
  {
    sparse_ ? sparse_matrix_.mtv(xx, yy) : dense_matrix_.mtv(xx, yy);
  }

  inline void add_to_entry(const size_t rr, const size_t cc, const ScalarType& value)
  {
    sparse_ ? sparse_matrix_.add_to_entry(rr, cc, value) : dense_matrix_.add_to_entry(rr, cc, value);
  }

  inline ScalarType get_entry(const size_t rr, const size_t cc) const
  {
    return sparse_ ? sparse_matrix_.get_entry(rr, cc) : dense_matrix_.get_entry(rr, cc);
  }

  inline void set_entry(const size_t rr, const size_t cc, const ScalarType value)
  {
    sparse_ ? sparse_matrix_.set_entry(rr, cc, value) : dense_matrix_.set_entry(rr, cc, value);
  }

  inline void clear_row(const size_t rr)
  {
    sparse_ ? sparse_matrix_.clear_row(rr) : dense_matrix_.clear_row(rr);
  }

  inline void clear_col(const size_t cc)
  {
    sparse_ ? sparse_matrix_.clear_col(cc) : dense_matrix_.clear_col(cc);
  }

  inline void unit_row(const size_t rr)
  {
    clear_row(rr);
    set_entry(rr, rr, ScalarType(1));
  }

  inline void unit_col(const size_t cc)
  {
    clear_col(cc);
    set_entry(cc, cc, ScalarType(1));
  }

  bool valid() const
  {
    return sparse_ ? sparse_matrix_.valid() : dense_matrix_.valid();
  }

  virtual size_t non_zeros() const override final
  {
    return sparse_ ? sparse_matrix_.non_zeros() : dense_matrix_.non_zeros();
  }

  virtual SparsityPatternDefault
  pattern(const bool prune = false,
          const EpsType eps = Common::FloatCmp::DefaultEpsilon<ScalarType>::value()) const override
  {
    return sparse_ ? sparse_matrix_.pattern(prune, eps) : dense_matrix_.pattern(prune, eps);
  } // ... pattern(...)

  /// \}

  template <class DuneDenseMatrixImp>
  void copy_to_densematrix(DuneDenseMatrixImp& ret) const
  {
    sparse_ ? sparse_matrix_.copy_to_densematrix(ret) : dense_matrix_.copy_to_densematrix(ret);
  }

  bool sparse() const
  {
    return sparse_;
  }

  SparseMatrixType& sparse_matrix()
  {
    return sparse_matrix_;
  }

  const SparseMatrixType& sparse_matrix() const
  {
    return sparse_matrix_;
  }

  DenseMatrixType& dense_matrix()
  {
    return dense_matrix_;
  }

  const DenseMatrixType& dense_matrix() const
  {
    return dense_matrix_;
  }

  template <class MatrixType>
  void rightmultiply(const MatrixType& other)
  {
    sparse_ ? sparse_matrix_.rightmultiply(other) : dense_matrix_.rightmultiply(other);
  } // void rightmultiply(...)

  void rightmultiply(const ThisType& other)
  {
    if (other.sparse())
      sparse_ ? sparse_matrix_.rightmultiply(other.sparse_matrix())
              : dense_matrix_.rightmultiply(other.sparse_matrix());
    else
      sparse_ ? sparse_matrix_.rightmultiply(other.dense_matrix()) : dense_matrix_.rightmultiply(other.dense_matrix());
  } // void rightmultiply(...)

  using MatrixInterfaceType::operator+;
  using MatrixInterfaceType::operator-;
  using MatrixInterfaceType::operator+=;
  using MatrixInterfaceType::operator-=;

  size_t num_rows_, num_cols_;
  bool sparse_;
  mutable SparseMatrixType sparse_matrix_;
  mutable DenseMatrixType dense_matrix_;
}; // class CommonSparseOrDenseMatrix<...>

template <class ScalarType = double>
using CommonSparseMatrixCsr = CommonSparseMatrix<ScalarType, SparseFormat::csr>;

template <class ScalarType = double>
using CommonSparseMatrixCsc = CommonSparseMatrix<ScalarType, SparseFormat::csc>;

template <class ScalarType = double>
using CommonSparseOrDenseMatrixCsr =
    CommonSparseOrDenseMatrix<CommonDenseMatrix<ScalarType>, CommonSparseMatrixCsr<ScalarType>>;

template <class ScalarType = double>
using CommonSparseOrDenseMatrixCsc =
    CommonSparseOrDenseMatrix<CommonDenseMatrix<ScalarType>, CommonSparseMatrixCsc<ScalarType>>;


} // namespace LA
namespace Common {


template <class T>
struct MatrixAbstraction<LA::CommonSparseMatrixCsr<T>>
    : public LA::internal::MatrixAbstractionBase<LA::CommonSparseMatrixCsr<T>>
{
};

template <class T>
struct MatrixAbstraction<LA::CommonSparseMatrixCsc<T>>
    : public LA::internal::MatrixAbstractionBase<LA::CommonSparseMatrixCsc<T>>
{
};

template <class DenseMatrixImp, class SparseMatrixImp>
struct MatrixAbstraction<LA::CommonSparseOrDenseMatrix<DenseMatrixImp, SparseMatrixImp>>
    : public LA::internal::MatrixAbstractionBase<LA::CommonSparseOrDenseMatrix<DenseMatrixImp, SparseMatrixImp>>
{
};


} // namespace Common
} // namespace XT
} // namespace Dune


// begin: this is what we need for the lib
#if DUNE_XT_WITH_PYTHON_BINDINGS


extern template class Dune::XT::LA::CommonSparseMatrix<double, SparseFormat::csr>;
extern template class Dune::XT::LA::CommonSparseMatrix<double, SparseFormat::csc>;


#endif // DUNE_XT_WITH_PYTHON_BINDINGS
// end: this is what we need for the lib


#endif // DUNE_XT_LA_CONTAINER_COMMON_MATRIX_SPARSE_HH
