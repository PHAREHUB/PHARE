#ifndef PHARE_PYTHON_PARTICLES_H
#define PHARE_PYTHON_PARTICLES_H

#include "pybind_def.h"
#include "kul/defs.hpp"

namespace PHARE::pydata
{
template<typename T, typename PyArray>
decltype(auto) make_ptr(PyArray const& py_array)
{
    auto ar_info = py_array.request();
    return kul::Pointers{static_cast<T*>(ar_info.ptr), static_cast<size_t>(ar_info.shape[0])};
}

template<size_t dim, typename PyArrayTuple>
decltype(auto) contiguous_view_from(PyArrayTuple const& py_particles)
{
    constexpr auto OwnedState = false;
    return core::ContiguousParticles<dim, OwnedState>{
        make_ptr<int32_t>(std::get<0>(py_particles)), // iCell
        make_ptr<float>(std::get<1>(py_particles)),   // delta
        make_ptr<double>(std::get<2>(py_particles)),  // weight
        make_ptr<double>(std::get<3>(py_particles)),  // charge
        make_ptr<double>(std::get<4>(py_particles))}; // v
}

template<size_t dim>
decltype(auto) make_py_array_tuple(size_t const size)
{
    return std::make_tuple(py_array_t<int32_t>(size * dim), // iCell
                           py_array_t<float>(size * dim),   // delta
                           py_array_t<double>(size),        // weight
                           py_array_t<double>(size),        // charge
                           py_array_t<double>(size * 3));   // v
}



template<size_t dim, typename PyArrayTuple>
void assert_particle_py_array_sizes(PyArrayTuple const& py_particles)
{
    auto shape_nd = [](auto const& ar, size_t dimdex = 0) { return ar.request().shape[dimdex]; };

    auto const& n_iCell  = shape_nd(std::get<0>(py_particles));
    auto const& n_delta  = shape_nd(std::get<1>(py_particles));
    auto const& n_weight = shape_nd(std::get<2>(py_particles));
    auto const& n_charge = shape_nd(std::get<3>(py_particles));

    assert(n_weight == n_charge);
    assert(n_delta == n_weight * dim);
    assert(n_iCell == n_weight * dim);

    auto const& v     = std::get<4>(py_particles);
    auto const& vDims = v.request().ndim;

    assert((vDims == 1 and shape_nd(v) == n_weight * 3)
           or (vDims == 2 and shape_nd(v) * shape_nd(v, 1) == n_weight * 3));
}

template<typename Splitter>
pyarray_particles_t split_pyarray_particles(pyarray_particles_crt const& py_particles)
{
    constexpr auto dim           = Splitter::dimension;
    constexpr auto nbRefinedPart = Splitter::nbRefinedPart;

    KUL_DEBUG_DO( // doesn't exist for release builds
        assert_particle_py_array_sizes<dim>(py_particles));

    auto particlesInView  = contiguous_view_from<dim>(py_particles);
    auto particlesOut     = make_py_array_tuple<dim>(particlesInView.size() * nbRefinedPart);
    auto particlesOutView = contiguous_view_from<dim>(particlesOut);

    Splitter splitter;
    for (size_t i = 0; i < particlesInView.size(); i++)
        splitter(particlesInView[i], particlesOutView, i * nbRefinedPart);

    return particlesOut;
}

} // namespace PHARE::pydata

#endif /*PHARE_PYTHON_PARTICLES_H*/
