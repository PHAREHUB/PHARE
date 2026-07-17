#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_DEF_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_DEF_HPP


#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/point/point.hpp"
#include "core/data/particles/particle.hpp"


#include <array>
#include <cstdint>

namespace PHARE::core
{


enum class LayoutMode : std::uint16_t {
    AoS = 0,
    AoSMapped, // 1
    AoSPC,     // 2
    AoSTS,     // 3
    AoSCMTS,   // 4
    SoA,       // 5
    SoAVX,     // 6
    SoATS,     // 7
    SoAVXTS,   // 8
    SoAPC,     // 9
    AoSPCTS,
};

bool constexpr is_tiled(LayoutMode lm)
{
    using enum LayoutMode;
    return any_in(lm, AoSTS, AoSCMTS, SoATS, SoAVXTS, AoSPCTS);
}

enum class StorageMode : std::uint16_t { ARRAY = 0, VECTOR, SPAN };

enum class ParticleType : std::uint16_t { Domain = 0, Ghost, PatchGhost, LevelGhost, All };


struct ParticleArrayOptions
{
    std::size_t dim          = 1;
    LayoutMode layout_mode   = LayoutMode::AoSMapped;
    StorageMode storage_mode = StorageMode::VECTOR;
    AllocatorMode alloc_mode = AllocatorMode::CPU;
    bool _const_             = 0; // sometimes needed

    auto constexpr with_layout(LayoutMode const lm) const
    {
        auto copy        = *this;
        copy.layout_mode = lm;
        return copy;
    }
    auto constexpr with_storage(StorageMode const sm) const
    {
        auto copy         = *this;
        copy.storage_mode = sm;
        return copy;
    }
    auto constexpr with_alloc(AllocatorMode const am) const
    {
        auto copy       = *this;
        copy.alloc_mode = am;
        return copy;
    }
};



template<std::size_t dim>
struct ParticleDefaults
{
    using Particle_t = Particle<dim>;
};


// carries where a particle came from when registering a move with move_check;
// per-cell layouts read icell (old AMR cell)
template<std::size_t dim>
struct ParticleTracker
{
    std::array<int, dim> icell{};
};

// tiled layouts additionally read tile_cell (old local tile cell)
template<std::size_t dim>
struct TiledParticleTracker
{
    std::array<int, dim> icell{};
    std::array<std::uint32_t, dim> tile_cell{};
};

// level ghost moves resolve differently: ghost cells are clamp-owned by border tiles,
// leaving the ghost box means deletion — distinct type so move_check can dispatch
template<std::size_t dim>
struct TiledLevelGhostParticleTracker
{
    std::array<int, dim> icell{};
    std::array<std::uint32_t, dim> tile_cell{};
};


// the tracker kind is a compile-time function of the layout and the particle type
template<LayoutMode layout_mode, ParticleType particle_type, std::size_t dim>
using ParticleTracker_t = std::conditional_t<
    is_tiled(layout_mode),
    std::conditional_t<particle_type == ParticleType::LevelGhost,
                       TiledLevelGhostParticleTracker<dim>, TiledParticleTracker<dim>>,
    ParticleTracker<dim>>;

template<LayoutMode layout_mode, ParticleType particle_type, std::size_t dim>
auto make_particle_tracker(auto&&... args) _PHARE_ALL_FN_
{
    return ParticleTracker_t<layout_mode, particle_type, dim>{args...};
}


template<typename R = std::uint32_t, std::size_t dim>
auto as_local_cell(std::array<int, dim> const& lower,
                   std::array<int, dim> const& icell) _PHARE_ALL_FN_
{
    return array_minus<R>(icell, lower);
}
template<std::size_t dim>
auto as_local_cell(Box<int, dim> const& box, std::array<int, dim> const& icell) _PHARE_ALL_FN_
{
    return as_local_cell(box.lower.toArray(), icell);
}
template<std::size_t dim>
auto as_local_cell(Box<int, dim> const& box, Point<int, dim> const& icell) _PHARE_ALL_FN_
{
    return as_local_cell(box.lower.toArray(), icell.toArray());
}

template<typename I0, typename I1> // support const vs non-const iterators
auto it_dist(I0&& i0, I1&& i1)
{
    auto const& begin = i0;
    auto const& pos   = i1;
    auto d            = std::distance(begin, pos);
    PHARE_ASSERT(d < 1e18); // should never happen // eg 2635249153387078728
    return d;
}



template<typename Box_t, typename RValue = std::uint32_t>
class LocalisedCellFlattener
{
public:
    static constexpr std::size_t dim = Box_t::dimension;

    LocalisedCellFlattener(Box_t const& b)
        : box{b}
        , shape{box.shape()}
    {
    }

    template<typename int_t>
    RValue operator()(std::array<int_t, dim> const icell) const
    {
        for (std::size_t i = 0; i < dim; ++i)
            icell[i] -= box.lower[i];
        if constexpr (dim == 2)
            return icell[1] + icell[0] * shape[1];
        if constexpr (dim == 3)
            return icell[2] + icell[1] * shape[2] + icell[0] * shape[1] * shape[2];
        return icell[0];
    }
    template<typename Particle>
    RValue operator()(Particle const& particle) const
    {
        return (*this)(particle.iCell);
    }

    Box_t const& box;

private:
    Point<int, dim> const& shape;
};




template<typename Particles_t, typename T, std::size_t D>
Particles_t make_particles(Box<T, D> const& box, std::size_t const ghost_cells)
{
    static_assert(Particles_t::storage_mode == StorageMode::VECTOR);

    using enum LayoutMode;
    if constexpr (any_in(Particles_t::layout_mode, AoSPC, SoAPC))
        return Particles_t{box, ghost_cells};
    else if constexpr (is_tiled(Particles_t::layout_mode))
        return Particles_t{box, ghost_cells};
    else if constexpr (any_in(Particles_t::layout_mode, AoSMapped))
        return Particles_t{grow(box, ghost_cells)};

    else
        return Particles_t{};
}

template<typename Particles_t, typename GridLayout_t>
Particles_t make_particles(GridLayout_t const& layout)
{
    return make_particles<Particles_t>(layout.AMRBox(), GridLayout_t::nbrParticleGhosts());
}



} // namespace PHARE::core


#endif /*PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_DEF_HPP*/
