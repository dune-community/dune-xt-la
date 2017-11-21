// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2014 - 2016)
//   Tobias Leibner  (2014 - 2015, 2017)

#include <dune/xt/common/test/main.hxx> // <- has to come first (includes the config.h)!

#include "eigensolver.hh"


struct EigenSolverForMatrixFrom2dEulerExample : public EigenSolverTestForMatricesWithRealEigenvaluesAndVectors
{
  using BaseType = EigenSolverTestForMatricesWithRealEigenvaluesAndVectors;
  using typename BaseType::MatrixType;
  using typename BaseType::ComplexMatrixType;
  using typename BaseType::RealMatrixType;
  using typename BaseType::EigenValuesType;
  using typename BaseType::RealEigenValuesType;

  /**
     From the literature, the known eigenvalues for the matrix below are
        "[-4.1419344967907001e-18 -4.1419344967907001e-18 0.7475656680523487 -0.7475656680523487]"
     and the known eigenvectors are
        "[1 0 2.6611618623465914 2.6611618623465914
          -0.0027456203140305249 3.9787864908411241 -0.0073065400681821043 -0.0073065400681821043
          -4.1419344967907001e-18 0 1.9893932454205621 -1.9893932454205621
          3.769215454408539e-06 -0.010924237014443618 3.7180152568215501 3.7180152568215501"

     Numpy (using lapack) gives the following eigenvalues
        "[7.47565668e-01 -7.47565668e-01 -4.38053098e-18 -4.13867822e-18]"
     and the following eigenvectors
        "[-5.33695006e-01  5.33695006e-01  9.99927103e-01 -9.99999986e-01;"
        "  1.46532385e-03 -1.46532385e-03  1.20742763e-02  1.64793482e-04;"
        " -3.98972064e-01 -3.98972064e-01 -4.40181522e-18  4.13836669e-18;"
        " -7.45646555e-01  7.45646555e-01 -3.69203189e-05  3.31675507e-06]"

     Lapack (this will be our reference below) gives the following eigenvalues
        (0.74756566805234859,0), (-0.7475656680523487,0), (-4.3805309824106515e-18,0), (-4.1386782167250046e-18,0)
     and the following eigenvectors
        "[-0.53369500623787047    0.53369500623787058    0.99992710258781181    -0.99999998641605381; "
        "  0.0014653238506238881 -0.0014653238506239037  0.012074276261282692    0.00016479348161139382; "
        " -0.39897206387455164   -0.39897206387455159   -4.4018152160265555e-18  4.1383666863327496e-18; "
        " -0.7456465552729743     0.74564655527297441   -3.6920318868549761e-05  3.3167550724756643e-06]");

     Sadly (bug report pending), eigen gives the following eigenvalues
        (0.7475656680523487,0), (-0.7475656680523487,0), (-4.1737887087823295e-18,6.7931406144174092e-20),
        (-4.1737887087823295e-18,-6.7931406144174092e-20)
     and the following eigenvectors
        (0.53369500623800015,0) (-0.53369500623800015,0) (-0.90295029334476462,0.42974388441822792)
        (-0.90295029334476462,-0.42974388441822792)
        (-0.0014653238506237894,0) (0.0014653238506237452,0) (5.6258735721227275e-06,0.00098055699157421558)
        (5.6258735721227275e-06,-0.00098055699157421558)
        (0.39897206387451295,0) (0.39897206387451295,0) (-2.6019023116072928e-17,-8.6866569199352384e-19)
        (-2.6019023116072928e-17,8.6866569199352384e-19)
        (0.74564655527290236,0) (-0.74564655527290247,0) (3.3879676871831931e-06,-4.312034485652854e-06)
        (3.3879676871831931e-06,4.312034485652854e-06)

   * \sa http://eigen.tuxfamily.org/bz/show_bug.cgi?id=1488
   */
  EigenSolverForMatrixFrom2dEulerExample()
  {
    matrix_ = XT::Common::from_string<MatrixType>(
        "[                      0                       0                       1                       0;"
        " -1.1372179493772346e-20 -4.1419344967907001e-18  -0.0027456203140305249                      -0;"
        "  1.5076861817634153e-06   0.0010982481256122097 -6.6270951948651204e-18     0.39999999999999991;"
        "  5.7868554526731795e-18  -4.548871797508937e-21      1.3971398393418408 -5.7987082955069798e-18]");
    expected_real_eigenvalues_ = XT::Common::from_string<RealEigenValuesType>(
        "[0.74756566805234859 -0.7475656680523487 -4.3805309824106515e-18 -4.1386782167250046e-18]");
    expected_real_eigenvectors_ = XT::Common::from_string<RealMatrixType>(
        "[-0.53369500623787047    0.53369500623787058    0.99992710258781181    -0.99999998641605381; "
        "  0.0014653238506238881 -0.0014653238506239037  0.012074276261282692    0.00016479348161139382; "
        " -0.39897206387455164   -0.39897206387455159   -4.4018152160265555e-18  4.1383666863327496e-18; "
        " -0.7456465552729743     0.74564655527297441   -3.6920318868549761e-05  3.3167550724756643e-06]");
    expected_eigenvalues_ = XT::Common::convert_to<EigenValuesType>(expected_real_eigenvalues_);
    expected_eigenvectors_ = XT::Common::convert_to<ComplexMatrixType>(expected_real_eigenvectors_);
    expected_max_ev_ = 0.74756566805234859;
    expected_min_ev_ = -0.7475656680523487;
    all_matrices_and_expected_eigenvalues_and_vectors_are_computed_ = true;
  }

  using BaseType::all_matrices_and_expected_eigenvalues_and_vectors_are_computed_;
  using BaseType::matrix_;
  using BaseType::expected_eigenvalues_;
  using BaseType::expected_eigenvectors_;
  using BaseType::expected_real_eigenvalues_;
  using BaseType::expected_max_ev_;
  using BaseType::expected_min_ev_;
  using BaseType::expected_real_eigenvectors_;
}; // struct EigenSolverForMatrixFrom2dEulerExample


TEST_F(EigenSolverForMatrixFrom2dEulerExample, exports_correct_types)
{
  exports_correct_types();
}

TEST_F(EigenSolverForMatrixFrom2dEulerExample, has_types_and_options)
{
  has_types_and_options();
}

TEST_F(EigenSolverForMatrixFrom2dEulerExample, throws_on_broken_matrix_construction)
{
  throws_on_broken_matrix_construction();
}

TEST_F(EigenSolverForMatrixFrom2dEulerExample, allows_broken_matrix_construction_when_checks_disabled)
{
  allows_broken_matrix_construction_when_checks_disabled();
}

TEST_F(EigenSolverForMatrixFrom2dEulerExample, throws_on_inconsistent_given_options)
{
  throws_on_inconsistent_given_options();
}

TEST_F(EigenSolverForMatrixFrom2dEulerExample, is_constructible)
{
  is_constructible();
}

TEST_F(EigenSolverForMatrixFrom2dEulerExample, gives_correct_eigenvalues)
{
  gives_correct_eigenvalues();
}

TEST_F(EigenSolverForMatrixFrom2dEulerExample, gives_correct_eigenvalues_in_correct_order)
{
  gives_correct_eigenvalues_in_correct_order();
}

TEST_F(EigenSolverForMatrixFrom2dEulerExample, gives_correct_real_eigenvalues)
{
  gives_correct_real_eigenvalues();
}

TEST_F(EigenSolverForMatrixFrom2dEulerExample, gives_correct_real_eigenvalues_in_correct_order)
{
  gives_correct_real_eigenvalues_in_correct_order();
}

TEST_F(EigenSolverForMatrixFrom2dEulerExample, gives_correct_max_eigenvalue)
{
  gives_correct_max_eigenvalue();
}

TEST_F(EigenSolverForMatrixFrom2dEulerExample, gives_correct_min_eigenvalue)
{
  gives_correct_min_eigenvalue();
}

TEST_F(EigenSolverForMatrixFrom2dEulerExample, gives_correct_eigenvectors_in_correct_order)
{
  gives_correct_eigenvectors_in_correct_order({{"eigen", /*we_expect_a_failure: */ "-1"}});
}

TEST_F(EigenSolverForMatrixFrom2dEulerExample, gives_correct_real_eigenvectors_in_correct_order)
{
  gives_correct_real_eigenvectors_in_correct_order({{"eigen", /*we_expect_a_failure: */ "-1"}});
}

TEST_F(EigenSolverForMatrixFrom2dEulerExample, gives_correct_eigendecomposition)
{
  gives_correct_eigendecomposition(1e-12);
}

TEST_F(EigenSolverForMatrixFrom2dEulerExample, gives_correct_real_eigendecomposition)
{
  gives_correct_real_eigendecomposition({{"lapack", "1e-12"}, {"eigen", /*we_expect_a_failure: */ "-1"}});
}