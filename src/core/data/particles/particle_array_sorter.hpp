#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SORTER_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SORTER_HPP

#include "core/def.hpp"
#include "core/data/vector.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/utilities/box/box.hpp"

#include <stdexcept>
#include <type_traits>

namespace PHARE::core
{

// --- public API ---

namespace detail
{
    // hides the per layout/alloc impl behind operator()/by_delta() - see implementation section
    // below
    template<auto layout_mde, auto alloc_mde, typename ParticleArray>
    struct ParticleSorter;
} // namespace detail

template<typename ParticleArray>
struct ParticleArraySorter
{
    using LM        = LayoutMode;
    using box_t     = Box<int, ParticleArray::dimension>;
    using sort_impl = detail::ParticleSorter<ParticleArray::layout_mode, ParticleArray::alloc_mode,
                                             ParticleArray>;

    bool constexpr static is_cell_sortable()
    {
        return !any_in(ParticleArray::layout_mode, LM::AoSPC, LM::AoSPCTS);
    }

    auto constexpr static cell_sortable = is_cell_sortable();

    void operator()()
    {
        if constexpr (cell_sortable)
            sort_impl{particles, box}();
    }

    void by_delta() { sort_impl{particles, box}.by_delta(); }

    ParticleArray& particles;
    box_t box;
};


struct Sortings // to be used in constexpr fashion
{
    bool by_delta = true;
};


template<Sortings S = Sortings{}, typename ParticleArray_t, typename Box_t>
auto& sort_particles(ParticleArray_t&& ps, Box_t const& box)
{
    using Particles = std::decay_t<ParticleArray_t>;

    if (ps.size() == 0)
        return ps;

    ParticleArraySorter<Particles> sorter{ps, box};
    sorter();
    if constexpr (S.by_delta)
        sorter.by_delta();

    if constexpr (Particles::layout_mode == LayoutMode::AoSPC)
        ps.reset_index_wrapper_map();

    return ps;
}


} // namespace PHARE::core


// --- implementations ---

namespace PHARE::core::detail
{

template<auto layout_mde, auto alloc_mde, typename ParticleArray>
struct ParticleSorter
{
    static_assert(all_are<LayoutMode>(layout_mde));
    static_assert(all_are<AllocatorMode>(alloc_mde));
    static_assert(dependent_false_v<ParticleArray>,
                  "ParticleSorter not implemented for this layout/alloc permutation");
};


using enum LayoutMode;
using enum AllocatorMode;

// shared by the flat (non-tiled, non-per-cell) layouts: AoS, AoSMapped, AoSPCTS
template<typename ParticleArray>
struct FlatParticleSorter
{
    using box_t = Box<int, ParticleArray::dimension>;

    FlatParticleSorter(ParticleArray& particles_, box_t const& domain)
        : particles{particles_}
        , domain_box{domain}
    {
    }

    ParticleArray& particles;
    box_t domain_box;
    CellFlattener<box_t> cell_flattener{domain_box};

    auto& operator()(std::int64_t const& l, std::int64_t const& r) // basically quicksort
    {
        std::sort(particles.begin() + l, particles.begin() + r,
                  [cf = cell_flattener](auto const& a, auto const& b) {
                      return cf(a.iCell()) < cf(b.iCell());
                  });
        return *this;
    }

    auto& operator()()
    {
        (*this)(0, particles.size());
        return *this;
    }

    auto static constexpr by_deltas()
    {
        return [](auto const& a, auto const& b) -> bool {
            return as_tuple(a.delta()) < as_tuple(b.delta());
        };
    }

    void by_deltas(std::uint64_t const& l, std::uint64_t const& r)
    {
        std::sort(particles.begin() + l, particles.begin() + r, by_deltas());
    }

    auto& by_delta()
    {
        // assumes already sorted by icell
        if (particles.size() == 0)
            return *this;

        auto const end = particles.end();
        auto beg       = particles.begin();
        auto lst       = particles.begin();

        auto const check = [&]() { return lst != end and lst.iCell() == beg.iCell(); };

        while (lst != end)
        {
            lst = beg + 1;
            while (check())
                ++lst;
            auto const s = it_dist(particles.begin(), beg);
            auto const e = it_dist(particles.begin(), lst);
            by_deltas(s, e);
            beg = lst;
        }

        return *this;
    }
};


template<typename ParticleArray>
struct ParticleSorter<AoS, CPU, ParticleArray> : FlatParticleSorter<ParticleArray>
{
    using FlatParticleSorter<ParticleArray>::FlatParticleSorter;
};

template<typename ParticleArray>
struct ParticleSorter<AoSMapped, CPU, ParticleArray> : FlatParticleSorter<ParticleArray>
{
    using FlatParticleSorter<ParticleArray>::FlatParticleSorter;
};

template<typename ParticleArray>
struct ParticleSorter<AoSPCTS, CPU, ParticleArray> : FlatParticleSorter<ParticleArray>
{
    using FlatParticleSorter<ParticleArray>::FlatParticleSorter;
};


// never actually invoked (ParticleArraySorter::is_cell_sortable() excludes AoSPC), kept for
// symmetry so forming a ParticleSorter<AoSPC, ...> fails at runtime rather than compile time
template<typename ParticleArray>
struct ParticleSorter<AoSPC, CPU, ParticleArray>
{
    using box_t = Box<int, ParticleArray::dimension>;

    auto& operator()(std::int64_t const& = 0, std::int64_t const& = 0)
    {
        throw std::runtime_error("no sort for AoSPC");
    }
    auto& operator()() { return (*this)(0, 0); }
    void by_deltas(std::uint64_t const&, std::uint64_t const&)
    {
        throw std::runtime_error("no sort for AoSPC");
    }
    auto& by_delta()
    {
        for (auto const& bix : particles.local_box())
            std::sort(particles(bix).begin(), particles(bix).end(),
                      [](auto const& a, auto const& b) -> bool {
                          return as_tuple(a.delta()) < as_tuple(b.delta());
                      });
        return *this;
    }

    ParticleArray& particles;
    box_t domain_box;
};


// SoA is internal-only plumbing (HDF5/restart write, python interop) and is never sorted -
// kept for symmetry, same reasoning as AoSPC above
template<typename ParticleArray>
struct ParticleSorter<SoA, CPU, ParticleArray>
{
    using box_t = Box<int, ParticleArray::dimension>;

    auto& operator()(std::int64_t const& = 0, std::int64_t const& = 0)
    {
        throw std::runtime_error("no sort for SoA");
    }
    auto& operator()() { return (*this)(0, 0); }
    void by_deltas(std::uint64_t const&, std::uint64_t const&)
    {
        throw std::runtime_error("no sort for SoA");
    }
    auto& by_delta() { throw std::runtime_error("no sort for SoA"); }

    ParticleArray& particles;
    box_t domain_box;
};


} // namespace PHARE::core::detail

#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SORTER_HPP */
