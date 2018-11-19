// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2018)

#ifndef DUNE_XT_LA_SOLVER_PBH
#define DUNE_XT_LA_SOLVER_PBH

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/operators.h>

#include <python/dune/xt/common/configuration.hh>
#include <dune/xt/common/python.hh>
#include <python/dune/xt/la/container.bindings.hh>

#include <dune/xt/la/container.hh>
#include <dune/xt/la/type_traits.hh>
#include <dune/xt/la/solver.hh>

//#define DBG_THROW(idx) DUNE_THROW(InvalidStateException, std::string("SEQ apply") + to_string(idx))
#define DBG_THROW(idx)

namespace Dune {
namespace XT {
namespace LA {


template <class M, class V = typename Container<typename M::ScalarType, M::vector_type>::VectorType>
typename std::enable_if<is_matrix<M>::value, void>::type bind_Solver(pybind11::module& m)
{
  using SeqComm = XT::SequentialCommunication;
  typedef Solver<M, SeqComm> C;
  using ParaComm = Dune::OwnerOverlapCopyCommunication<unsigned long, int>;
  typedef Solver<M, ParaComm> ParaC;

  namespace py = pybind11;
  using namespace pybind11::literals;

  const auto ClassName = Common::to_camel_case(bindings::container_name<M>::value() + "_solver");
  py::class_<C> c(m, ClassName.c_str(), ClassName.c_str());

  c.def_static("types", &C::types);
  c.def_static("options", &C::options);

  c.def(py::init<M>());

  c.def("apply",
        [](const C& self, const V& rhs, V& solution) {
          DBG_THROW(1);
          self.apply(rhs, solution);
        },
        "rhs"_a,
        "solution"_a);
  c.def("apply",
        [](const C& self, const V& rhs, V& solution, const std::string& type) {
          DBG_THROW(2);
          self.apply(rhs, solution, type);
        },
        "rhs"_a,
        "solution"_a,
        "type"_a);
  c.def("apply",
        [](const C& self, const V& rhs, V& solution, const Common::Configuration& options) {
          DBG_THROW(3);
          self.apply(rhs, solution, options);
        },
        "rhs"_a,
        "solution"_a,
        "options"_a);


  const auto ParaClassName = Common::to_camel_case(bindings::container_name<M>::value() + "_solver_parallel");
  py::class_<ParaC> parac(m, ParaClassName.c_str(), ParaClassName.c_str());

  parac.def_static("types", &C::types);
  parac.def_static("options", &C::options);

  parac.def(py::init<M, const ParaComm&>());

  parac.def(
      "apply", [](const ParaC& self, const V& rhs, V& solution) { self.apply(rhs, solution); }, "rhs"_a, "solution"_a);
  parac.def(
      "apply",
      [](const ParaC& self, const V& rhs, V& solution, const std::string& type) { self.apply(rhs, solution, type); },
      "rhs"_a,
      "solution"_a,
      "type"_a);
  parac.def("apply",
            [](const ParaC& self, const V& rhs, V& solution, const Common::Configuration& options) {
              self.apply(rhs, solution, options);
            },
            "rhs"_a,
            "solution"_a,
            "options"_a);

  m.def("make_solver",
        [](const M& matrix) {
          DBG_THROW(4);
          return C(matrix);
        },
        pybind11::keep_alive<0, 1>());
  m.def(
      "make_solver", [](const M& matrix, const SeqComm& /*comm*/) { return C(matrix); }, pybind11::keep_alive<0, 1>());

  m.def("make_solver",
        [](const M& matrix, const ParaComm& dof_communicator) { return ParaC(matrix, dof_communicator); },
        pybind11::keep_alive<0, 1>());

} // ... bind_Solver(...)


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_SOLVER_PBH
