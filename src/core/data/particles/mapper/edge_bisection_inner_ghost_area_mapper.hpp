#ifndef PHARE_CORE_DATA_PARTICLES_EDGE_BISECTION_INNER_GHOST_MAPPER_HPP
#define PHARE_CORE_DATA_PARTICLES_EDGE_BISECTION_INNER_GHOST_MAPPER_HPP

#include <cstddef>

#include "core/logger.hpp"
#include "core/vector.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/range/range.hpp"

namespace PHARE::core
{
template<typename ParticleArray, std::size_t ghost_particle_width_>
class EdgeBisectionMapper
{
    static constexpr auto dim                  = ParticleArray::dimension;
    static constexpr auto ghost_particle_width = ghost_particle_width_;

    using Ranges = std::vector<std::pair<std::size_t*, std::size_t*>>;

public:
    template<typename Box_t>
    auto map(ParticleArray const& /*particles*/, Box_t const& /*box*/);
};

template<typename Box_t, typename Icell>
auto flat_cell(Box_t const& box, Icell const& icell)
{
    if constexpr (Box_t::dimension == 3)
        return icell[2] + icell[1] * box.upper[2] + icell[0] * box.upper[1] * box.upper[0];

    return 1;
}


template<typename Ps, typename Box_t, typename Icell>
auto find_left_of(Ps const& particles, Box_t const& box, Icell const& icell, std::size_t idx)
{
    // if(){}

    // if(particles)
    PHARE_LOG_LINE_STR("find_left_of");

    return 1;
}

template<typename Ps, typename Box_t, typename Icell>
auto find_right_of(Ps const& particles, Box_t const& box, Icell const& icell, std::size_t idx)
{
    PHARE_LOG_LINE_STR("find_right_of");

    auto find_cell = [&]() {
        auto v = icell;
        v[0] += 1;
        return v;
    }();

    auto flat_cell0 = flat_cell(box, icell);
    auto flat_cell1 = flat_cell(box, find_cell);

    while (true)
    {
        if (flat_cell(box, particles.iCell(idx)) < flat_cell1)
        {
            auto half = idx / 2;
            idx       = +()
        }
    }

    return idx;
}

template<std::size_t ghost_particle_width_, typename Ps, typename Box_t, typename Icell>
auto bisect_left(Ps const& ps, Box_t const& box, Icell const& icell, std::size_t lo, std::size_t up)
{
    auto const left_wall_icell = [&]() {
        auto v = icell;
        v[0]   = box.lower[0];
        return v;
    }();
    auto const left_ghost_max = [&]() {
        auto v = left_wall_icell;
        v[0] += ghost_particle_width_;
        return v;
    }();

    auto left_of_lower  = find_left_of(particles, box, left_wall_icell, idx);
    auto right_of_upper = find_right_of(particles, box, left_ghost_max, idx);

    return 1;
}

template<std::size_t ghost_particle_width_, typename ParticleArray, typename Box_t, typename Icell>
auto bisect_right(ParticleArray const& particles, Box_t const& box, Icell const& icell)
{
}



template<typename ParticleArray, std::size_t ghost_particle_width_>
template<typename Box_t>
auto EdgeBisectionMapper<ParticleArray, ghost_particle_width_>::map(ParticleArray const& particles,
                                                                    Box_t const& box)
{
    static_assert(dim == 3); // for now

    Ranges ranges;

    std::size_t mid_index    = particles.size() / 2;
    auto const& mid_particle = particles.iCell(mid_index);

    bisect_left<ghost_particle_width_>(particles, box, mid_particle, mid_index);
    bisect_right<ghost_particle_width_>(particles, box, mid_particle);

    for (int z = 0; z <= box.upper[2]; ++z)
        for (int y = 0; z <= box.upper[1]; ++y)
        {
            //
        }
}

} // namespace PHARE::core
#endif /* PHARE_CORE_DATA_PARTICLES_EDGE_BISECTION_INNER_GHOST_MAPPER_HPP */
