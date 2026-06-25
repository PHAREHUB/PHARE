#ifndef PHARE_PYTHON_PYBIND_DEF_HPP
#define PHARE_PYTHON_PYBIND_DEF_HPP

#include <tuple>
#include <cassert>
#include <cstdint>
#include <stdexcept>

#include "core/utilities/span.hpp"
#include "core/utilities/types.hpp"

#include "pybind11/stl.h"
#include "pybind11/numpy.h"


namespace PHARE::pydata
{
template<typename T>
using py_array_t = pybind11::array_t<T, pybind11::array::c_style | pybind11::array::forcecast>;


using pyarray_particles_t = std::tuple<py_array_t<int32_t>, py_array_t<double>, py_array_t<double>,
                                       py_array_t<double>, py_array_t<double>>;

using pyarray_particles_crt
    = std::tuple<py_array_t<int32_t> const&, py_array_t<double> const&, py_array_t<double> const&,
                 py_array_t<double> const&, py_array_t<double> const&>;

template<typename PyArrayInfo>
std::size_t ndSize(PyArrayInfo const& ar_info)
{
    assert(ar_info.ndim >= 1 && ar_info.ndim <= 3);

    return std::accumulate(ar_info.shape.begin(), ar_info.shape.end(), 1,
                           std::multiplies<std::size_t>());
}


template<std::size_t dim>
auto shape(pybind11::buffer_info const& info)
{
    if (info.ndim != dim)
        throw std::runtime_error("bad dim");
    return core::for_N_make_array<dim>([&](auto i) -> std::size_t { return info.shape[i]; });
}

template<std::size_t dim, typename T>
auto shape(py_array_t<T> const& ar)
{
    return shape<dim>(ar.request());
}


template<std::size_t dim>
auto strides(pybind11::buffer_info const& info)
{
    if (info.ndim != dim)
        throw std::runtime_error("bad dim");
    return core::for_N_make_array<dim>([&](auto i) -> std::size_t { return info.strides[i]; });
}

template<std::size_t dim, typename T>
auto strides(py_array_t<T> const& ar)
{
    return strides<dim>(ar.request());
}

template<typename T, typename Size, std::size_t dim>
std::array<std::size_t, dim> strides_from(std::array<Size, dim> const& shape)
{
    std::size_t constexpr data_size = sizeof(T);

    if constexpr (dim == 1)
        return {data_size};

    if constexpr (dim == 2)
        return {data_size * shape[1], data_size};

    if constexpr (dim == 3)
        return {data_size * shape[1] * shape[2], data_size * shape[1], data_size};
}

template<typename T>
class __attribute__((visibility("hidden"))) PyArrayWrapper : public core::Span<T>
{
public:
    PyArrayWrapper(py_array_t<T> const& array)
        : core::Span<T>{static_cast<T*>(array.request().ptr), pydata::ndSize(array.request())}
        , _array{array}
    {
        assert(_array.request().ptr);
        assert(_array.request().ptr == array.request().ptr); // assert no copy
    }

protected:
    py_array_t<T> _array;
};

template<typename T>
std::shared_ptr<core::Span<T>> makePyArrayWrapper(py_array_t<T> const& array)
{
    return std::make_shared<PyArrayWrapper<T>>(array);
}

template<typename T>
core::Span<T> makeSpan(py_array_t<T> const& py_array)
{
    auto ar_info = py_array.request();
    assert(ar_info.ptr);
    return {static_cast<T*>(ar_info.ptr), ndSize(ar_info)};
}


template<typename Field>
auto field_as_memory_view(Field& field)
{
    using value_type = Field::value_type;
    return pybind11::memoryview::from_buffer(field.data(), field.shape(),
                                             strides_from<value_type>(field.shape()));
}


} // namespace PHARE::pydata

#endif /*PHARE_PYTHON_PYBIND_DEF_H*/
