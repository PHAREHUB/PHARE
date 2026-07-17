#ifndef PHARE_CORE_DATA_PARTICLES_SORTING_PARTICLES_SORTING_HPP
#define PHARE_CORE_DATA_PARTICLES_SORTING_PARTICLES_SORTING_HPP

#include "core/def.hpp"
#include "core/vector.hpp"
#include "core/utilities/box/box.hpp"

#include "core/data/particles/arrays/particle_array_soa_thrust.hpp"


namespace PHARE::core::detail
{
template<auto alloc_mode, std::size_t impl, typename ParticleArray>
class ParticleSorter;
} // namespace PHARE::core::detail

#include "particles_sorting_cpu.hpp"

#if PHARE_HAVE_MKN_GPU
#include "particles_sorting_gpu_mkn.hpp"
#else // no impl (yet)
namespace PHARE::core::detail
{
template<typename ParticleArray>
class ParticleSorter<AllocatorMode::GPU_UNIFIED, /*impl = */ 0, ParticleArray>;
// template<typename ParticleArray>
// class ParticleSorter<alloc_mode::GPU, /*impl = */ 0, ParticleArray>;
} // namespace PHARE::core::detail
#endif



namespace PHARE::core
{



template<typename ParticleArray, std::size_t impl = 0>
struct ParticleArraySorter
{
    using LM        = LayoutMode;
    using box_t     = Box<int, ParticleArray::dimension>;
    using sort_impl = detail::ParticleSorter<ParticleArray::alloc_mode, impl, ParticleArray>;

    bool constexpr static is_sortable()
    {
        return true; //! any_in(ParticleArray::layout_mode, AoSTS, SoATS);
    }

    bool constexpr static is_cell_sortable()
    {
        return !any_in(ParticleArray::layout_mode, LM::AoSPC, LM::SoAPC);
    }

    auto constexpr static sortable      = is_sortable();
    auto constexpr static cell_sortable = is_cell_sortable();

    void operator()()
    {
        if constexpr (sortable and cell_sortable)
            sort_impl{particles, box}();
    }

    void by_delta()
    {
        if constexpr (sortable)
            sort_impl{particles, box}.by_delta();
    }


    ParticleArray& particles;
    box_t box;
};

} // namespace PHARE::core




#endif /*PHARE_CORE_DATA_PARTICLES_SORTING_PARTICLES_SORTING_HPP*/
