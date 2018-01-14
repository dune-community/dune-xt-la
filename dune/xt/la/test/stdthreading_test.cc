// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-common developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#include <config.h>

#include <memory>
#include <thread>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#endif

#include <dune/xt/common/string.hh>
#include <dune/xt/la/container/pattern.hh>
#include <dune/xt/la/container/istl.hh>

using namespace Dune;


int main(int argc, char** argv)
{
  if (argc < 5)
    return EXIT_FAILURE;

#if HAVE_DUNE_FEM
  Fem::MPIManager::initialize(argc, argv);
#else
  MPIHelper::instance(argc, argv);
#endif

  const auto N = XT::Common::from_string<size_t>(argv[1]);
  const auto S = XT::Common::from_string<size_t>(argv[2]);
  const auto M = XT::Common::from_string<size_t>(argv[3]);
  const auto W = XT::Common::from_string<size_t>(argv[4]);

  std::cout << "computing " << N << "x" << N << " unit matrix with " << S << " entries per row ... " << std::flush;
  Timer timer;

  XT::LA::SparsityPatternDefault pattern(N);
  for (size_t ii = 0; ii < N; ++ii)
    pattern.insert(ii, ii);
  for (size_t jj = 1; jj < S; ++jj) {
    for (size_t ii = jj; ii < N; ++ii)
      pattern.insert(ii, ii - jj);
    for (size_t ii = 0; ii < N - jj; ++ii)
      pattern.insert(ii, ii + jj);
  }
  pattern.sort();

  XT::LA::IstlRowMajorSparseMatrix<double> mat(N, N, pattern);
  for (size_t ii = 0; ii < N; ++ii)
    mat.unit_row(ii);

  std::cout << "done (took " << timer.elapsed() << "s)" << std::endl;
  std::cout << "preparing " << M << " input vectors ... " << std::flush;
  timer.reset();

  std::vector<std::shared_ptr<XT::LA::IstlDenseVector<double>>> Us(M);
  std::vector<std::shared_ptr<XT::LA::IstlDenseVector<double>>> Vs(M);
  for (size_t ii = 0; ii < M; ++ii) {
    Us[ii] = std::make_shared<XT::LA::IstlDenseVector<double>>(N, ii);
    Vs[ii] = std::make_shared<XT::LA::IstlDenseVector<double>>(N, 0);
  }

  std::cout << "done (took " << timer.elapsed() << "s)" << std::endl;
  std::cout << "doing mv with " << W << " threads ... " << std::flush;
  timer.reset();

  std::vector<size_t> partition_starts(W);
  for (size_t w = 0; w < W; ++w)
    partition_starts[w] = w * (M / W);
  std::vector<size_t> partition_ends(W);
  for (size_t w = 0; w < W - 1; ++w)
    partition_ends[w] = partition_starts[w + 1];
  partition_ends[W - 1] = M;
  std::vector<std::shared_ptr<std::thread>> threads(W, nullptr);
  for (size_t w = 0; w < W; ++w)
    threads[w] = std::make_shared<std::thread>(
        [&](const size_t p_s, const size_t p_e) {
          for (size_t ii = p_s; ii < p_e; ++ii)
            mat.mv(*Us[ii], *Vs[ii]);
        },
        partition_starts[w],
        partition_ends[w]);
  for (size_t w = 0; w < W; ++w)
    threads[w]->join();

  std::cout << "done (took " << timer.elapsed() << "s)" << std::endl;

  for (size_t ii = 0; ii < M; ++ii)
    if (Us[ii]->sup_norm() != Vs[ii]->sup_norm())
      DUNE_THROW(InvalidStateException, ii);
} // ... main(...)
