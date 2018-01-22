# This file is part of the dune-xt-la project:
#   https://github.com/dune-community/dune-xt-la
# Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2017)
#   Rene Milk       (2017 - 2018)
#   Tobias Leibner  (2017)

from dune.xt.codegen import typeid_to_typedef_name as safe_name, have_eigen

matrix = ['EigenDenseMatrix<double>', 'FieldMatrix<double, 2, 2>']
field = ['double', 'double']
complex_matrix = ['EigenDenseMatrix<std::complex<double>>', 'FieldMatrix<std::complex<double>, 2, 2>']
real_matrix = ['EigenDenseMatrix<double>', 'FieldMatrix<double, 2, 2>']


def _ok(ft):
    if 'Eigen' in ft[0]:
        return have_eigen(cache)
    return True

testtypes = [(safe_name('_'.join(ft)), *ft) for ft in zip(matrix, field, complex_matrix, real_matrix) if _ok(ft)]