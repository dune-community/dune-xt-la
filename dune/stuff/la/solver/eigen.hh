// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff/
// Copyright Holders: Felix Albrecht, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_LA_SOLVER_EIGEN_HH
#define DUNE_STUFF_LA_SOLVER_EIGEN_HH

#if HAVE_EIGEN

#include <string>
#include <vector>
#include <algorithm>
#include <sstream>

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#if HAVE_UMFPACK
#include <Eigen/UmfPackSupport>
#endif
//#include <Eigen/SPQRSupport>
//#include <Eigen/CholmodSupport>
#if HAVE_SUPERLU
#include <Eigen/SuperLUSupport>
#endif

#if HAVE_FASP
extern "C" {
#include "fasp_functs.h"
}
#endif // HAVE_FASP

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/configtree.hh>
#include <dune/stuff/la/container/eigen.hh>

#include "../solver.hh"

namespace Dune {
namespace Stuff {
namespace LA {


template <class S>
class Solver<EigenDenseMatrix<S>> : protected SolverUtils
{
public:
  typedef EigenDenseMatrix<S> MatrixType;

  Solver(const MatrixType& matrix)
    : matrix_(matrix)
  {
  }

  static std::vector<std::string> options()
  {
    return {"qr.colpivhouseholder",
            "qr.fullpivhouseholder",
            "qr.householder",
            "lu.fullpiv",
            "llt",
            "ldlt",
            "lu.partialpiv"};
  }

  static Common::ConfigTree options(const std::string& type)
  {
    SolverUtils::check_given(type, options());
    return Common::ConfigTree("type", type);
  } // ... options(...)

  size_t apply(const EigenDenseVector<S>& rhs, EigenDenseVector<S>& solution) const
  {
    return redirect_to_appropriate_apply(rhs, solution);
  }

  size_t apply(const EigenDenseVector<S>& rhs, EigenMappedDenseVector<S>& solution) const
  {
    return redirect_to_appropriate_apply(rhs, solution);
  }

  size_t apply(const EigenMappedDenseVector<S>& rhs, EigenDenseVector<S>& solution) const
  {
    return redirect_to_appropriate_apply(rhs, solution);
  }

  size_t apply(const EigenMappedDenseVector<S>& rhs, EigenMappedDenseVector<S>& solution) const
  {
    return redirect_to_appropriate_apply(rhs, solution);
  }

  size_t apply(const EigenDenseVector<S>& rhs, EigenDenseVector<S>& solution, const std::string& type) const
  {
    return redirect_to_appropriate_apply(rhs, solution, type);
  }

  size_t apply(const EigenDenseVector<S>& rhs, EigenMappedDenseVector<S>& solution, const std::string& type) const
  {
    return redirect_to_appropriate_apply(rhs, solution, type);
  }

  size_t apply(const EigenMappedDenseVector<S>& rhs, EigenDenseVector<S>& solution, const std::string& type) const
  {
    return redirect_to_appropriate_apply(rhs, solution, type);
  }

  size_t apply(const EigenMappedDenseVector<S>& rhs, EigenMappedDenseVector<S>& solution, const std::string& type) const
  {
    return redirect_to_appropriate_apply(rhs, solution, type);
  }

  size_t apply(const EigenDenseVector<S>& rhs, EigenDenseVector<S>& solution, const Common::ConfigTree& opts) const
  {
    return redirect_to_appropriate_apply(rhs, solution, opts);
  }

  size_t apply(const EigenDenseVector<S>& rhs, EigenMappedDenseVector<S>& solution,
               const Common::ConfigTree& opts) const
  {
    return redirect_to_appropriate_apply(rhs, solution, opts);
  }

  size_t apply(const EigenMappedDenseVector<S>& rhs, EigenDenseVector<S>& solution,
               const Common::ConfigTree& opts) const
  {
    return redirect_to_appropriate_apply(rhs, solution, opts);
  }

  size_t apply(const EigenMappedDenseVector<S>& rhs, EigenMappedDenseVector<S>& solution,
               const Common::ConfigTree& opts) const
  {
    return redirect_to_appropriate_apply(rhs, solution, opts);
  }

private:
  template <class RhsType, class SolutionType>
  size_t redirect_to_appropriate_apply(const RhsType& rhs, SolutionType& solution) const
  {
    return apply(rhs, solution, options()[0]);
  } // ... redirect_to_appropriate_apply(...)

  template <class RhsType, class SolutionType>
  size_t redirect_to_appropriate_apply(const RhsType& rhs, SolutionType& solution, const Common::ConfigTree& opts) const
  {
    if (!opts.has_key("type"))
      DUNE_THROW_COLORFULLY(Exception::configuration_error,
                            "Given options (see below) need to have at least the key 'type' set!\n\n" << opts);
    const auto type = opts.get<std::string>("type");
    return apply(rhs, solution, type);
  } // ... redirect_to_appropriate_apply(...)

  template <class RhsType, class SolutionType>
  size_t redirect_to_appropriate_apply(const RhsType& rhs, SolutionType& solution, const std::string& type) const
  {
    this->check_given(type, options());
    // solve
    if (type == "qr.colpivhouseholder") {
      solution.backend() = matrix_.backend().colPivHouseholderQr().solve(rhs.backend());
    } else if (type == "qr.fullpivhouseholder")
      solution.backend() = matrix_.backend().fullPivHouseholderQr().solve(rhs.backend());
    else if (type == "qr.householder")
      solution.backend() = matrix_.backend().householderQr().solve(rhs.backend());
    else if (type == "lu.fullpiv")
      solution.backend() = matrix_.backend().fullPivLu().solve(rhs.backend());
    else if (type == "llt")
      solution.backend() = matrix_.backend().llt().solve(rhs.backend());
    else if (type == "ldlt")
      solution.backend() = matrix_.backend().ldlt().solve(rhs.backend());
    else if (type == "lu.partialpiv")
      solution.backend() = matrix_.backend().partialPivLu().solve(rhs.backend());
    else
      DUNE_THROW_COLORFULLY(Exception::internal_error,
                            "Given type '" << type << "' is not supported, although it was reported by options()!");
// check
#ifndef NDEBUG
    if (!rhs.backend().isApprox(matrix_.backend() * solution.backend()))
      return 4;
#endif // NDEBUG
    return 0;
  } // ... redirect_to_appropriate_apply(...)

  const MatrixType& matrix_;
}; // class Solver


/**
 *  \note lu.sparse will copy the matrix to column major
 *  \note qr.sparse will copy the matrix to column major
 *  \note ldlt.simplicial will copy the matrix to column major
 *  \note llt.simplicial will copy the matrix to column major
 *  \note llt.simplicial will copy the matrix to column major
 *  \note spqr will copy the matrix to column major
 */
template <class S>
class Solver<EigenRowMajorSparseMatrix<S>> : protected SolverUtils
{
  typedef ::Eigen::SparseMatrix<S, ::Eigen::ColMajor> ColMajorBackendType;

public:
  typedef EigenRowMajorSparseMatrix<S> MatrixType;


  Solver(const MatrixType& matrix)
    : matrix_(matrix)
  {
  }

  static std::vector<std::string> options()
  {
    return
    {
      "bicgstab.ilut", "bicgstab.diagonal", "bicgstab.identity", "lu.sparse", "qr.sparse", "ldlt.simplicial",
          "llt.simplicial" /*, "spqr", "llt.cholmodsupernodal"*/
#if HAVE_UMFPACK
          ,
          "lu.umfpack"
#endif
#if HAVE_SUPERLU
          ,
          "superlu"
#endif
          ,
          "cg.diagonal.lower", "cg.diagonal.upper", "cg.identity.lower", "cg.identity.upper"
    };
  }

  static Common::ConfigTree options(const std::string& type)
  {
    // check
    SolverUtils::check_given(type, options());
    // direct solvers
    if (type == "lu.sparse" || type == "qr.sparse" || type == "ldlt.simplicial" || type == "llt.simplicial"
        || type == "lu.umfpack"
        || type == "spqr"
        || type == "llt.cholmodsupernodal"
        || type == "superlu")
      return Common::ConfigTree("type", type);
    // iterative solvers
    Common::ConfigTree iterative_options({"max_iter", "precision"}, {"10000", "1e-10"});
    if (type == "bicgstab.ilut") {
      iterative_options.add("preconditioner.fill_factor", "10");
      iterative_options.add("preconditioner.drop_tol", "1e-4");
    }
    iterative_options.add("type", type);
    return iterative_options;
  } // ... options(...)

  size_t apply(const EigenDenseVector<S>& rhs, EigenDenseVector<S>& solution) const
  {
    return redirect_to_appropriate_apply(rhs, solution);
  }

  size_t apply(const EigenDenseVector<S>& rhs, EigenMappedDenseVector<S>& solution) const
  {
    return redirect_to_appropriate_apply(rhs, solution);
  }

  size_t apply(const EigenMappedDenseVector<S>& rhs, EigenDenseVector<S>& solution) const
  {
    return redirect_to_appropriate_apply(rhs, solution);
  }

  size_t apply(const EigenMappedDenseVector<S>& rhs, EigenMappedDenseVector<S>& solution) const
  {
    return redirect_to_appropriate_apply(rhs, solution);
  }

  size_t apply(const EigenDenseVector<S>& rhs, EigenDenseVector<S>& solution, const std::string& type) const
  {
    return redirect_to_appropriate_apply(rhs, solution, type);
  }

  size_t apply(const EigenDenseVector<S>& rhs, EigenMappedDenseVector<S>& solution, const std::string& type) const
  {
    return redirect_to_appropriate_apply(rhs, solution, type);
  }

  size_t apply(const EigenMappedDenseVector<S>& rhs, EigenDenseVector<S>& solution, const std::string& type) const
  {
    return redirect_to_appropriate_apply(rhs, solution, type);
  }

  size_t apply(const EigenMappedDenseVector<S>& rhs, EigenMappedDenseVector<S>& solution, const std::string& type) const
  {
    return redirect_to_appropriate_apply(rhs, solution, type);
  }

  size_t apply(const EigenDenseVector<S>& rhs, EigenDenseVector<S>& solution, const Common::ConfigTree& opts) const
  {
    return redirect_to_appropriate_apply(rhs, solution, opts);
  }

  size_t apply(const EigenDenseVector<S>& rhs, EigenMappedDenseVector<S>& solution,
               const Common::ConfigTree& opts) const
  {
    return redirect_to_appropriate_apply(rhs, solution, opts);
  }

  size_t apply(const EigenMappedDenseVector<S>& rhs, EigenDenseVector<S>& solution,
               const Common::ConfigTree& opts) const
  {
    return redirect_to_appropriate_apply(rhs, solution, opts);
  }

  size_t apply(const EigenMappedDenseVector<S>& rhs, EigenMappedDenseVector<S>& solution,
               const Common::ConfigTree& opts) const
  {
    return redirect_to_appropriate_apply(rhs, solution, opts);
  }

private:
  template <class RhsType, class SolutionType>
  size_t redirect_to_appropriate_apply(const RhsType& rhs, SolutionType& solution) const
  {
    return redirect_to_appropriate_apply(rhs, solution, options()[0]);
  }

  template <class RhsType, class SolutionType>
  size_t redirect_to_appropriate_apply(const RhsType& rhs, SolutionType& solution, const Common::ConfigTree& opts) const
  {
    if (!opts.has_key("type"))
      DUNE_THROW_COLORFULLY(Exception::configuration_error,
                            "Given options (see below) need to have at least the key 'type' set!\n\n" << opts);
    const auto type                       = opts.get<std::string>("type");
    const Common::ConfigTree default_opts = options(type);
    if (type == "cg.diagonal.lower") {
      // create solver
      typedef ::Eigen::ConjugateGradient<typename MatrixType::BackendType,
                                         ::Eigen::Lower,
                                         ::Eigen::DiagonalPreconditioner<S>> SolverType;
      SolverType solver(matrix_.backend());
      solver.setMaxIterations(opts.get("max_iter", default_opts.get<std::size_t>("max_iter")));
      solver.setTolerance(opts.get("precision", default_opts.get<S>("precision")));
      solution.backend() = solver.solve(rhs.backend());
      if (solver.info() != ::Eigen::Success)
        return solver.info();
    } else if (type == "cg.diagonal.upper") {
      typedef ::Eigen::ConjugateGradient<typename MatrixType::BackendType,
                                         ::Eigen::Upper,
                                         ::Eigen::DiagonalPreconditioner<double>> SolverType;
      SolverType solver(matrix_.backend());
      solver.setMaxIterations(opts.get("max_iter", default_opts.get<std::size_t>("max_iter")));
      solver.setTolerance(opts.get("precision", default_opts.get<S>("precision")));
      solution.backend() = solver.solve(rhs.backend());
      if (solver.info() != ::Eigen::Success)
        return solver.info();
    } else if (type == "cg.identity.lower") {
      typedef ::Eigen::ConjugateGradient<typename MatrixType::BackendType,
                                         ::Eigen::Lower,
                                         ::Eigen::IdentityPreconditioner> SolverType;
      SolverType solver(matrix_.backend());
      solver.setMaxIterations(opts.get("max_iter", default_opts.get<std::size_t>("max_iter")));
      solver.setTolerance(opts.get("precision", default_opts.get<S>("precision")));
      solution.backend() = solver.solve(rhs.backend());
      if (solver.info() != ::Eigen::Success)
        return solver.info();
    } else if (type == "cg.identity.upper") {
      typedef ::Eigen::ConjugateGradient<typename MatrixType::BackendType,
                                         ::Eigen::Lower,
                                         ::Eigen::IdentityPreconditioner> SolverType;
      SolverType solver(matrix_.backend());
      solver.setMaxIterations(opts.get("max_iter", default_opts.get<std::size_t>("max_iter")));
      solver.setTolerance(opts.get("precision", default_opts.get<S>("precision")));
      solution.backend() = solver.solve(rhs.backend());
      if (solver.info() != ::Eigen::Success)
        return solver.info();
    } else if (type == "bicgstab.ilut") {
      typedef ::Eigen::BiCGSTAB<typename MatrixType::BackendType, ::Eigen::IncompleteLUT<S>> SolverType;
      SolverType solver(matrix_.backend());
      solver.setMaxIterations(opts.get("max_iter", default_opts.get<std::size_t>("max_iter")));
      solver.setTolerance(opts.get("precision", default_opts.get<S>("precision")));
      solver.preconditioner().setDroptol(
          opts.get("preconditioner.drop_tol", default_opts.get<S>("preconditioner.drop_tol")));
      solver.preconditioner().setFillfactor(
          opts.get("preconditioner.fill_factor", default_opts.get<size_t>("preconditioner.fill_factor")));
      solution.backend() = solver.solve(rhs.backend());
      if (solver.info() != ::Eigen::Success)
        return solver.info();
    } else if (type == "bicgstab.diagonal") {
      typedef ::Eigen::BiCGSTAB<typename MatrixType::BackendType, ::Eigen::DiagonalPreconditioner<S>> SolverType;
      SolverType solver(matrix_.backend());
      solver.setMaxIterations(opts.get("max_iter", default_opts.get<std::size_t>("max_iter")));
      solver.setTolerance(opts.get("precision", default_opts.get<S>("precision")));
      solution.backend() = solver.solve(rhs.backend());
      if (solver.info() != ::Eigen::Success)
        return solver.info();
    } else if (type == "bicgstab.identity") {
      typedef ::Eigen::BiCGSTAB<typename MatrixType::BackendType, ::Eigen::IdentityPreconditioner> SolverType;
      SolverType solver(matrix_.backend());
      solver.setMaxIterations(opts.get("max_iter", default_opts.get<std::size_t>("max_iter")));
      solver.setTolerance(opts.get("precision", default_opts.get<S>("precision")));
      solution.backend() = solver.solve(rhs.backend());
      if (solver.info() != ::Eigen::Success)
        return solver.info();
    } else if (type == "lu.sparse") {
      ColMajorBackendType colmajor_copy(matrix_.backend());
      colmajor_copy.makeCompressed();
      typedef ::Eigen::SparseLU<ColMajorBackendType> SolverType;
      SolverType solver;
      solver.analyzePattern(colmajor_copy);
      solver.factorize(colmajor_copy);
      solution.backend() = solver.solve(rhs.backend());
      if (solver.info() != ::Eigen::Success)
        return solver.info();
    } else if (type == "qr.sparse") {
      ColMajorBackendType colmajor_copy(matrix_.backend());
      colmajor_copy.makeCompressed();
      typedef ::Eigen::SparseQR<ColMajorBackendType, ::Eigen::COLAMDOrdering<int>> SolverType;
      SolverType solver;
      solver.analyzePattern(colmajor_copy);
      solver.factorize(colmajor_copy);
      solution.backend() = solver.solve(rhs.backend());
      if (solver.info() != ::Eigen::Success)
        return solver.info();
    } else if (type == "ldlt.simplicial") {
      ColMajorBackendType colmajor_copy(matrix_.backend());
      colmajor_copy.makeCompressed();
      typedef ::Eigen::SimplicialLDLT<ColMajorBackendType> SolverType;
      SolverType solver;
      solver.analyzePattern(colmajor_copy);
      solver.factorize(colmajor_copy);
      solution.backend() = solver.solve(rhs.backend());
      if (solver.info() != ::Eigen::Success)
        return solver.info();
    } else if (type == "llt.simplicial") {
      ColMajorBackendType colmajor_copy(matrix_.backend());
      colmajor_copy.makeCompressed();
      typedef ::Eigen::SimplicialLLT<ColMajorBackendType> SolverType;
      SolverType solver;
      solver.analyzePattern(colmajor_copy);
      solver.factorize(colmajor_copy);
      solution.backend() = solver.solve(rhs.backend());
      if (solver.info() != ::Eigen::Success)
        return solver.info();
#if HAVE_UMFPACK
    } else if (type == "lu.umfpack") {
      typedef ::Eigen::UmfPackLU<typename MatrixType::BackendType> SolverType;
      SolverType solver;
      solver.analyzePattern(matrix_.backend());
      solver.factorize(matrix_.backend());
      solution.backend() = solver.solve(rhs.backend());
      if (solver.info() != ::Eigen::Success)
        return solver.info();
#endif // HAVE_UMFPACK
//    } else if (type == "spqr") {
//      ColMajorBackendType colmajor_copy(matrix_.backend());
//      colmajor_copy.makeCompressed();
//      typedef ::Eigen::SPQR< ColMajorBackendType > SolverType;
//      SolverType solver;
//      solver.analyzePattern(colmajor_copy);
//      solver.factorize(colmajor_copy);
//      solution.backend() = solver.solve(rhs.backend());
//      if (solver.info() != ::Eigen::Success)
//        return solver.info();
//    } else if (type == "cholmodsupernodalllt") {
//      typedef ::Eigen::CholmodSupernodalLLT< typename MatrixType::BackendType > SolverType;
//      SolverType solver;
//      solver.analyzePattern(matrix_.backend());
//      solver.factorize(matrix_.backend());
//      solution.backend() = solver.solve(rhs.backend());
//      if (solver.info() != ::Eigen::Success)
//        return solver.info();
#if HAVE_SUPERLU
    } else if (type == "superlu") {
      typedef ::Eigen::SuperLU<typename MatrixType::BackendType> SolverType;
      SolverType solver;
      solver.analyzePattern(matrix_.backend());
      solver.factorize(matrix_.backend());
      solution.backend() = solver.solve(rhs.backend());
      if (solver.info() != ::Eigen::Success)
        return solver.info();
#endif // HAVE_SUPERLU
    } else
      DUNE_THROW_COLORFULLY(Exception::internal_error,
                            "Given type '" << type << "' is not supported, although it was reported by options()!");
// check
#ifndef NDEBUG
    if (!rhs.backend().isApprox(matrix_.backend() * solution.backend()))
      return 4;
#endif // NDEBUG
    return 0;
  } // ... redirect_to_appropriate_apply(...)

  template <class RhsType, class SolutionType>
  size_t redirect_to_appropriate_apply(const RhsType& rhs, SolutionType& solution, const std::string& type) const
  {
    return redirect_to_appropriate_apply(rhs, solution, options(type));
  }

  const MatrixType& matrix_;
}; // class Solver


} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // HAVE_EIGEN

#endif // DUNE_STUFF_LA_SOLVER_EIGEN_HH
