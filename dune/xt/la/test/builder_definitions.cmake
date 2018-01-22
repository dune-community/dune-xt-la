# This file is part of the dune-xt-la project:
#   https://github.com/dune-community/dune-xt-la
# Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Rene Milk      (2017 - 2018)
#   Tobias Leibner (2018)

set(DXT_BIN_COUNT "13" CACHE STRING "number of bins for test targets" )
add_custom_target(test_binaries_builder_0 DEPENDS test_solver)
set_tests_properties(test_solver PROPERTIES LABELS "builder_0")
add_custom_target(test_binaries_builder_1 DEPENDS headercheck__dune_xt_la_container_common_vector_dense.hh headercheck__dune_xt_la_container_vector-interface.hh headercheck__dune_xt_la_eigen-solver_internal_lapacke.hh test_eigensolver_for_real_matrix_with_real_evs)
set_tests_properties(test_eigensolver_for_real_matrix_with_real_evs PROPERTIES LABELS "builder_1")
add_custom_target(test_binaries_builder_2 DEPENDS headercheck__dune_xt_la_container_io.hh headercheck__dune_xt_la_container_matrix-interface.hh test_eigensolver_for_real_matrix_with_real_evs_from_2d_euler_equations)
set_tests_properties(test_eigensolver_for_real_matrix_with_real_evs_from_2d_euler_equations PROPERTIES LABELS "builder_2")
add_custom_target(test_binaries_builder_3 DEPENDS headercheck__dune_xt_la_algorithms_cholesky.hh headercheck__dune_xt_la_matrix-inverter_internal_base.hh test_eigensolver_for_real_matrix_with_complex_evs)
set_tests_properties(test_eigensolver_for_real_matrix_with_complex_evs PROPERTIES LABELS "builder_3")
add_custom_target(test_binaries_builder_4 DEPENDS headercheck__dune_xt_la_container_conversion.hh headercheck__dune_xt_la_test_container.hh test_eigensolver_for_real_matrix_with_distinct_real_evs)
set_tests_properties(test_eigensolver_for_real_matrix_with_distinct_real_evs PROPERTIES LABELS "builder_4")
add_custom_target(test_binaries_builder_5 DEPENDS headercheck__dune_xt_la_container_common_matrix.hh headercheck__dune_xt_la_container_eigen.hh headercheck__dune_xt_la_solver_istl_amg.hh test_container)
set_tests_properties(test_container PROPERTIES LABELS "builder_5")
add_custom_target(test_binaries_builder_6 DEPENDS headercheck__dune_xt_la_algorithms_solve_sym_tridiag_posdef.hh headercheck__dune_xt_la_container.hh headercheck__dune_xt_la_container_common_vector_sparse.hh headercheck__dune_xt_la_container_eigen_base.hh headercheck__dune_xt_la_container_pattern.hh headercheck__dune_xt_la_exceptions.hh headercheck__dune_xt_la_type_traits.hh test_container_vector)
set_tests_properties(test_container_vector PROPERTIES LABELS "builder_6")
add_custom_target(test_binaries_builder_7 DEPENDS headercheck__dune_xt_la_container.bindings.hh headercheck__dune_xt_la_container_container-interface.hh headercheck__dune_xt_la_container_eigen_sparse.hh headercheck__dune_xt_la_eigen-solver_internal_numpy.hh headercheck__dune_xt_la_solver.hh)
add_custom_target(test_binaries_builder_8 DEPENDS headercheck__dune_xt_la_container_common_matrix_sparse.hh headercheck__dune_xt_la_container_eigen_dense.hh headercheck__dune_xt_la_container_eye-matrix.hh test_container_matrix test_empty)
set_tests_properties(test_container_matrix PROPERTIES LABELS "builder_8")
set_tests_properties(test_empty PROPERTIES LABELS "builder_8")
add_custom_target(test_binaries_builder_9 DEPENDS headercheck__dune_xt_la_algorithms.hh headercheck__dune_xt_la_algorithms_qr.hh headercheck__dune_xt_la_algorithms_triangular_solves.hh headercheck__dune_xt_la_matrix-inverter.hh)
add_custom_target(test_binaries_builder_10 DEPENDS headercheck__dune_xt_la_container_common_matrix_dense.hh headercheck__dune_xt_la_container_interfaces.hh headercheck__dune_xt_la_eigen-solver_eigen.hh headercheck__dune_xt_la_eigen-solver_fmatrix.hh headercheck__dune_xt_la_matrix-inverter_eigen.hh headercheck__dune_xt_la_matrix-inverter_internal_eigen.hh headercheck__dune_xt_la_solver_istl.hh)
add_custom_target(test_binaries_builder_11 DEPENDS headercheck__dune_xt_la_container_common_vector.hh headercheck__dune_xt_la_container_istl.hh headercheck__dune_xt_la_container_unit_matrices.hh headercheck__dune_xt_la_container_vector-interface-internal.hh headercheck__dune_xt_la_eigen-solver_internal_shifted-qr.hh headercheck__dune_xt_la_solver_common.hh headercheck__dune_xt_la_solver_fasp.hh headercheck__dune_xt_la_test_eigensolver.hh)
add_custom_target(test_binaries_builder_12 DEPENDS headercheck__dune_xt_la_container_common.hh headercheck__dune_xt_la_eigen-solver.hh headercheck__dune_xt_la_eigen-solver_internal_base.hh headercheck__dune_xt_la_eigen-solver_internal_eigen.hh headercheck__dune_xt_la_solver_eigen.hh)