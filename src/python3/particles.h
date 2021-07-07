#ifndef PHARE_PYTHON_PARTICLES_H
#define PHARE_PYTHON_PARTICLES_H

#include <cassert>
#include <cstddef>
#include <tuple>
#include "amr/data/particles/refine/particles_data_split.h"
#include "core/data/particles/particle_packer.h"
#include "core/data/particles/particle.h"
#include "core/utilities/types.h"

#include "python3/pybind_def.h"

namespace PHARE::pydata
{
template<std::size_t dim, typename PyArrayTuple>
core::ContiguousParticlesView<dim> contiguousViewFrom(PyArrayTuple const& py_particles)
{
    return {makeSpan<int>(std::get<0>(py_particles)),     // iCell
            makeSpan<double>(std::get<1>(py_particles)),  // delta
            makeSpan<double>(std::get<2>(py_particles)),  // weight
            makeSpan<double>(std::get<3>(py_particles)),  // charge
            makeSpan<double>(std::get<4>(py_particles))}; // v
}

template<std::size_t dim>
pyarray_particles_t makePyArrayTuple(std::size_t const size)
{
    return std::make_tuple(py_array_t<int>(size * dim),    // iCell
                           py_array_t<double>(size * dim), // delta
                           py_array_t<double>(size),       // weight
                           py_array_t<double>(size),       // charge
                           py_array_t<double>(size * 3));  // v
}



template<std::size_t dim, typename PyArrayTuple>
void assertParticlePyArraySizes(PyArrayTuple const& py_particles)
{
    auto shape_nd = [](auto const& ar, std::size_t dimdex = 0) -> std::size_t {
        return ar.request().shape[dimdex];
    };

    auto const& n_iCell  = ndSize(std::get<0>(py_particles).request());
    auto const& n_delta  = ndSize(std::get<1>(py_particles).request());
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
pyarray_particles_t splitPyArrayParticles(pyarray_particles_crt const& py_particles)
{
    constexpr auto dim           = Splitter::dimension;
    constexpr auto interp_order  = Splitter::interp_order;
    constexpr auto nbRefinedPart = Splitter::nbRefinedPart;

    PHARE_DEBUG_DO(assertParticlePyArraySizes<dim>(py_particles));

    auto particlesInView  = contiguousViewFrom<dim>(py_particles);
    auto particlesOut     = makePyArrayTuple<dim>(particlesInView.size() * nbRefinedPart);
    auto particlesOutView = contiguousViewFrom<dim>(particlesOut);

    Splitter splitter;

    for (std::size_t i = 0; i < particlesInView.size(); i++)
        splitter(amr::toFineGrid<interp_order>(std::copy(particlesInView[i])), particlesOutView,
                 i * nbRefinedPart);

    return particlesOut;
}

} // namespace PHARE::pydata

#endif /*PHARE_PYTHON_PARTICLES_H*/
