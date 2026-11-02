#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SORTER_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SORTER_HPP

#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/sorting/particles_sorting.hpp"
#include <type_traits>

namespace PHARE::core
{

struct Sortings // to be used in constexpr fashion
{
    bool by_delta = true;
};


template<Sortings S = Sortings{}, typename ParticleArray_t, typename Box_t>
auto& sort_particles(ParticleArray_t&& ps, Box_t const& box)
{
    using Particles                         = std::decay_t<ParticleArray_t>;
    std::string_view constexpr static FN_ID = "sort,";
    auto constexpr function_id              = join_string_views_v<FN_ID, Particles::type_id>;
    PHARE_LOG_SCOPE(1, function_id);

    if (ps.size() == 0)
        return ps;

    ParticleArraySorter<Particles> sorter{ps, box};
    sorter();
    if constexpr (S.by_delta)
        sorter.by_delta();

    using enum LayoutMode;
    if constexpr (any_in(Particles::layout_mode, AoSPC, SoAPC))
        ps.reset_index_wrapper_map();

    return ps;
}



} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SORTER_HPP */
