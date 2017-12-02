// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)

#include "config.h"

#if HAVE_DUNE_PYBINDXI

#include <string>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#endif

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/common/bindings.hh>

#include "container/container-interface.pbh"
#include "container/vector-interface.pbh"
#include "container/pattern.pbh"
#include "container/matrix-interface.pbh"
#include "solver.pbh"

#include "container.hh"


PYBIND11_PLUGIN(_la)
{
  namespace py = pybind11;
  using namespace pybind11::literals;
  namespace LA = Dune::XT::LA;

  py::module m("_la", "dune-xt-la");
  DUNE_XT_COMMON_BINDINGS_INITIALIZE(m, "dune.xt.la");

  LA::bind_Backends(m);

//  auto common_dense_vector_double = LA::bind_Vector<LA::CommonDenseVector<double>>(m);
#if HAVE_DUNE_ISTL
  auto istl_dense_vector_double = LA::bind_Vector<LA::IstlDenseVector<double>>(m);
#endif
  //#if HAVE_EIGEN
  //  auto eigen_dense_vector_double = LA::bind_Vector<LA::EigenDenseVector<double>>(m);
  //#endif

  LA::bind_SparsityPatternDefault(m);

#define BIND_MATRIX(C, s, c) auto c = LA::bind_Matrix<C, s>(m);

//  BIND_MATRIX(LA::CommonDenseMatrix<double>, false, common_dense_matrix_double);
//  BIND_MATRIX(LA::CommonSparseMatrix<double>, true, common_sparse_matrix_double);
#if HAVE_DUNE_ISTL
  BIND_MATRIX(LA::IstlRowMajorSparseMatrix<double>, true, istl_row_major_sparse_matrix_double);
#endif
//#if HAVE_EIGEN
//  BIND_MATRIX(LA::EigenDenseMatrix<double>, false, eigen_dense_matrix_double);
//  BIND_MATRIX(LA::EigenRowMajorSparseMatrix<double>, true, eigen_row_major_sparse_matrix_double);
//#endif
#undef BIND_MATRIX
//  LA::addbind_Matrix_Vector_interaction(common_dense_matrix_double, common_dense_vector_double);
//  LA::addbind_Matrix_Vector_interaction(common_sparse_matrix_double, common_dense_vector_double);
#if HAVE_DUNE_ISTL
  LA::addbind_Matrix_Vector_interaction(istl_row_major_sparse_matrix_double, istl_dense_vector_double);
#endif
//#if HAVE_EIGEN
//  LA::addbind_Matrix_Vector_interaction(eigen_dense_matrix_double, eigen_dense_vector_double);
//  LA::addbind_Matrix_Vector_interaction(eigen_row_major_sparse_matrix_double, eigen_dense_vector_double);
//#endif

//  LA::bind_Solver<LA::CommonDenseMatrix<double>>(m);
//  LA::bind_Solver<LA::CommonSparseMatrix<double>>(m);
#if HAVE_DUNE_ISTL
  LA::bind_Solver<LA::IstlRowMajorSparseMatrix<double>>(m);
#endif
  //#if HAVE_EIGEN
  //  LA::bind_Solver<LA::EigenDenseMatrix<double>>(m);
  //  LA::bind_Solver<LA::EigenRowMajorSparseMatrix<double>>(m);
  //#endif

  return m.ptr();
} // PYBIND11_PLUGIN(la)

#endif // HAVE_DUNE_PYBINDXI
