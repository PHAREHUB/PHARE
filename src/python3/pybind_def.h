#ifndef PHARE_PYTHON_PYBIND_DEF_H
#define PHARE_PYTHON_PYBIND_DEF_H

#include <tuple>
#include <cstdint>
#include <stdexcept>

#include "pybind11/stl.h"
#include "pybind11/numpy.h"

namespace PHARE::pydata
{
template<typename T>
using py_array_t = pybind11::array_t<T, pybind11::array::c_style | pybind11::array::forcecast>;


using pyarray_particles_t = std::tuple<py_array_t<int32_t>, py_array_t<float>, py_array_t<double>,
                                       py_array_t<double>, py_array_t<double>>;

using pyarray_particles_crt
    = std::tuple<py_array_t<int32_t> const&, py_array_t<float> const&, py_array_t<double> const&,
                 py_array_t<double> const&, py_array_t<double> const&>;



template<typename T>
core::Span<T, int> to_span(py_array_t<T> const& py_array)
{
    py::buffer_info info = py_array.request();
    if (!info.ptr)
        throw std::runtime_error("to_span received an array with an invalid internal ptr");
    assert(info.ndim == 1 or info.ndim == 2);
    int size = info.ndim == 1 ? info.shape[0] : info.shape[0] * info.shape[1];
    return {static_cast<T*>(info.ptr), size};
}


} // namespace PHARE::pydata

#endif /*PHARE_PYTHON_PYBIND_DEF_H*/
