// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2018)

#ifndef DUNE_XT_LA_CONTAINER_INTERFACE_PBH
#define DUNE_XT_LA_CONTAINER_INTERFACE_PBH

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/operators.h>

#include <dune/xt/la/type_traits.hh>

#include <dune/xt/la/container/container-interface.hh>

namespace Dune {
namespace XT {
namespace LA {


pybind11::enum_<Backends> bind_Backends(pybind11::module& m)
{
  namespace py = pybind11;

  py::enum_<Backends> c(m, "Backends");
  c.value("common_dense", Backends::common_dense);
  c.value("common_sparse", Backends::common_sparse);
  c.value("istl_sparse", Backends::istl_sparse);
  c.value("eigen_dense", Backends::eigen_dense);
  c.value("eigen_sparse", Backends::eigen_sparse);
  c.value("none", Backends::none);

  m.attr("default_backend") = py::cast(default_backend);
  m.attr("default_sparse_backend") = py::cast(default_sparse_backend);
  m.attr("default_dense_backend") = py::cast(default_dense_backend);

  return c;
} // ... bind_Backends(...)


template <class C>
typename std::enable_if<is_container<C>::value, void>::type addbind_ContainerInterface(pybind11::class_<C>& c)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  typedef typename C::ScalarType S;

  c.def("copy",
        [](C& self, const bool deep) {
          if (deep)
            return self.copy();
          else
            return C(self);
        },
        "deep"_a = false);
  c.def("scal", [](C& self, const S& alpha) { self.scal(alpha); }, "alpha"_a);
  c.def("axpy", [](C& self, const S& alpha, const C& xx) { self.axpy(alpha, xx); });
  c.def("has_equal_shape", [](const C& self, const C& other) { return self.has_equal_shape(other); }, "other"_a);
  c.def(py::self *= S());
  c.def(py::self * S());
  c.def("__sub__", [](const C& self, const C& other) {
    auto ret = self.copy();
    ret.axpy(-1, other);
    return ret;
  });

} // ... addbind_ContainerInterface(...)


template <class C>
typename std::enable_if<provides_backend<C>::value, void>::type addbind_ProvidesBackend(pybind11::class_<C>& c)
{
  namespace py = pybind11;

  c.def_property_readonly_static("backend_type", [](py::object /*self*/) { return C::backend_type; });
}

template <class C>
typename std::enable_if<!provides_backend<C>::value, void>::type addbind_ProvidesBackend(pybind11::class_<C>& /*c*/)
{}


/**
 * \brief Allows the resulting container to be convertible into a NumPy array as in `np.array(c, copy = False)`.
 */
template <class C>
typename std::enable_if<provides_data_access<C>::value, pybind11::class_<C>>::type
bind_ProvidesDataAccess(pybind11::module& m, const std::string& class_id, const std::string& help_id)
{
  namespace py = pybind11;
  typedef typename C::DataType D;

  py::class_<C> c(m, class_id.c_str(), help_id.c_str(), py::buffer_protocol());

  c.def_buffer([](C& vec) -> py::buffer_info {
    return py::buffer_info(
        vec.data(), sizeof(D), py::format_descriptor<D>::format(), 1, {vec.data_size()}, {sizeof(D)});
  });

  return c;
}

template <class C>
typename std::enable_if<!provides_data_access<C>::value, pybind11::class_<C>>::type
bind_ProvidesDataAccess(pybind11::module& m, const std::string& class_id, const std::string& help_id)
{
  namespace py = pybind11;
  return py::class_<C>(m, class_id.c_str(), help_id.c_str());
}


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_CONTAINER_INTERFACE_PBH
