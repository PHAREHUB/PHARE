#ifndef PHARE_PYTHON_PYBIND_DEF_H
#define PHARE_PYTHON_PYBIND_DEF_H

#include "phare/include.h"

#include "pybind11/stl.h"
#include "pybind11/numpy.h"
#include "pybind11/chrono.h"
#include "pybind11/complex.h"
#include "pybind11/functional.h"

namespace PHARE::pydata
{
template<typename T>
using py_array_t = pybind11::array_t<T, pybind11::array::c_style | pybind11::array::forcecast>;


using pyarray_particles_t = std::tuple<py_array_t<int32_t>, py_array_t<float>, py_array_t<double>,
                                       py_array_t<double>, py_array_t<double>>;

using pyarray_particles_crt
    = std::tuple<py_array_t<int32_t> const&, py_array_t<float> const&, py_array_t<double> const&,
                 py_array_t<double> const&, py_array_t<double> const&>;


} // namespace PHARE::pydata

#endif /*PHARE_PYTHON_PYBIND_DEF_H*/
