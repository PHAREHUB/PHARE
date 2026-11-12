// to #include "core/data/particles/particle_array_exporter.hpp"

#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_EXPORTER
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_EXPORTER

#include "core/def.hpp"
#include "core/def/phare_config.hpp" // IWYU pragma: keep

#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/box/box_spans.hpp"
#include "core/utilities/range/range.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_array_appender.hpp"
#include "core/data/particles/particle_array_partitioner.hpp"

#include <cstddef>
#include <iterator>
#include <algorithm>
#include <stdexcept>

namespace PHARE::core
{

// --- public API ---

// hides the per layout/alloc impl behind its methods - see implementation section below
template<auto layout_mde, auto alloc_mde>
struct ParticlesExporter;


template<typename Src, typename Box_t>
void delete_particles_in(Src& src, Box_t const& box)
{
    using Exporter = ParticlesExporter<Src::layout_mode, Src::alloc_mode>;

    Exporter{}.delete_particles_in(src, box);
}

template<typename Src, typename Box_t>
void delete_particles_not_in(Src& src, Box_t const& box)
{
    using Exporter = ParticlesExporter<Src::layout_mode, Src::alloc_mode>;

    Exporter{}.delete_particles_not_in(src, box);
}


template<typename Dst, typename Src, typename Box_t>
void move_in_domain(Dst& dst, Src& src, Box_t const& domain_box)
{
    static_assert(Dst::layout_mode == Src::layout_mode);
    using Exporter = ParticlesExporter<Src::layout_mode, Src::alloc_mode>;
    Exporter{}.move_in_domain(dst, src, domain_box);
}

template<typename Dst, typename Src, typename T, std::size_t dim> //
void move_in_ghost_layer(Dst& dst, Src& src, Box<T, dim> const& domain_box,
                         Box<T, dim> const& ghost_box)
{
    static_assert(Dst::layout_mode == Src::layout_mode);
    using Exporter = ParticlesExporter<Src::layout_mode, Src::alloc_mode>;
    Exporter{}.move_in_ghost_layer(dst, src, domain_box, ghost_box);
}

template<typename Dst, typename Src, typename T, std::size_t dim, typename Boxes>
void move_in_ghost_layer(Dst& dst, Src& src, Box<T, dim> const& domain_box,
                         Boxes const& ghost_boxes)
{
    static_assert(Dst::layout_mode == Src::layout_mode);
    using Exporter = ParticlesExporter<Src::layout_mode, Src::alloc_mode>;
    Exporter{}.move_in_ghost_layer(dst, src, domain_box, ghost_boxes);
}


// --- implementations ---

template<auto layout_mde, auto alloc_mde>
struct ParticlesExporter
{
    static_assert(all_are<LayoutMode>(layout_mde));
    static_assert(all_are<AllocatorMode>(alloc_mde));

    auto constexpr static layout_mode = layout_mde;
    auto constexpr static alloc_mode  = alloc_mde;

    template<typename Src, typename Dst, typename Box_t>
    void move_particles(Src& src, Dst& dst, Box_t const& box, std::size_t const growby = 0);


    template<typename Src, std::size_t dim>
    void delete_particles_in(Src& src, Box<int, dim> const& box);

    template<typename Src, typename Boxes>
    void delete_particles_in(Src& src, Boxes const& boxes);

    template<typename Src, std::size_t dim>
    void delete_particles_not_in(Src& src, Box<int, dim> const& box);

    template<typename Src, typename Boxes>
    void delete_particles_not_in(Src& src, Boxes const& boxes);

    template<typename Dst, typename Src, std::size_t dim>
    void move_in_domain(Dst& dst, Src& src, Box<int, dim> const& domain_box);

    template<typename Dst, typename Src, std::size_t dim>
    void move_in_ghost_layer(Dst& dst, Src& src, Box<int, dim> const& domain_box,
                             Box<int, dim> const& ghost_box);

    template<typename Dst, typename Src, std::size_t dim, typename Boxes>
    void move_in_ghost_layer(Dst& dst, Src& src, Box<int, dim> const& domain_box,
                             Boxes const& ghost_boxes);

    // template<typename Src, typename Dst, typename Box_t, typename Refiner, typename Transformer>
    // void operator()(Src const&, Dst&, Box_t const&, Refiner, Transformer);
};


using enum LayoutMode;
using enum AllocatorMode;


template<>
template<typename Src, std::size_t dim>
void ParticlesExporter<AoS, CPU>::delete_particles_not_in(Src& src, Box<int, dim> const& box)
{
    auto const kept = std::partition(src.begin(), src.end(),
                                     [&](auto const& p) { return isIn(p.iCell(), box); });
    src.erase(kept, src.end());
}

template<>
template<typename Src, typename Boxes>
void ParticlesExporter<AoS, CPU>::delete_particles_not_in(Src& src, Boxes const& boxes)
{
    auto const kept = std::partition(src.begin(), src.end(),
                                     [&](auto const& p) { return isIn(p.iCell(), boxes); });
    src.erase(kept, src.end());
}

template<>
template<typename Src, typename Boxes>
void ParticlesExporter<AoSMapped, CPU>::delete_particles_not_in(Src& src, Boxes const& boxes)
{
    auto range      = makeIndexRange(src);
    auto const kept = src.partition(range, [&](auto const& cell) { return isIn(cell, boxes); });
    src.erase(makeRange(src, kept.iend(), src.size()));
}

template<>
template<typename Src, std::size_t dim>
void ParticlesExporter<AoSMapped, CPU>::delete_particles_not_in(Src& src, Box<int, dim> const& box)
{
    auto range      = makeIndexRange(src);
    auto const kept = src.partition(range, [&](auto const& cell) { return isIn(cell, box); });
    src.erase(makeRange(src, kept.iend(), src.size()));
}

template<>
template<typename Src, typename Boxes>
void ParticlesExporter<AoS, CPU>::delete_particles_in(Src& src, Boxes const& boxes)
{
    auto const kept = std::partition(src.begin(), src.end(),
                                     [&](auto const& p) { return !isIn(p.iCell(), boxes); });
    src.erase(kept, src.end());
}

template<>
template<typename Dst, typename Src, std::size_t dim>
void ParticlesExporter<AoS, CPU>::move_in_domain(Dst& dst, Src& src,
                                                 Box<int, dim> const& domain_box)
{
    auto const moved = std::partition(src.begin(), src.end(),
                                      [&](auto const& p) { return !isIn(p.iCell(), domain_box); });
    std::copy(moved, src.end(), std::back_inserter(dst));
    src.erase(moved, src.end());
}

template<>
template<typename Dst, typename Src, std::size_t dim>
void ParticlesExporter<AoS, CPU>::move_in_ghost_layer(Dst& dst, Src& src,
                                                      Box<int, dim> const& domain_box,
                                                      Box<int, dim> const& ghost_box)
{
    auto const ghost_layer_boxes = ghost_box.remove(domain_box);
    auto const moved             = std::partition(
        src.begin(), src.end(), [&](auto const& p) { return !isIn(p.iCell(), ghost_layer_boxes); });
    std::copy(moved, src.end(), std::back_inserter(dst));
    src.erase(moved, src.end());
}

template<>
template<typename Dst, typename Src, std::size_t dim, typename Boxes>
void ParticlesExporter<AoS, CPU>::move_in_ghost_layer(Dst& dst, Src& src,
                                                      Box<int, dim> const& domain_box,
                                                      Boxes const& ghost_boxes)
{
    for (auto const& gb : ghost_boxes)
        this->move_in_ghost_layer(dst, src, domain_box, gb);
}

template<>
template<typename Src, typename Boxes>
void ParticlesExporter<AoSMapped, CPU>::delete_particles_in(Src& src, Boxes const& boxes)
{
    auto range      = makeIndexRange(src);
    auto const kept = src.partition(range, [&](auto const& cell) { return !isIn(cell, boxes); });
    src.erase(makeRange(src, kept.iend(), src.size()));
}

// AoSCMTS: one flat cell-mapped array per tile, see MappedTileSetVector::erase_if
template<>
template<typename Src, std::size_t dim>
void ParticlesExporter<AoSCMTS, CPU>::delete_particles_not_in(Src& src, Box<int, dim> const& box)
{
    src.erase_if([&](auto const& p) { return !isIn(p.iCell(), box); });
}

template<>
template<typename Src, typename Boxes>
void ParticlesExporter<AoSCMTS, CPU>::delete_particles_not_in(Src& src, Boxes const& boxes)
{
    src.erase_if([&](auto const& p) { return !isIn(p.iCell(), boxes); });
}

template<>
template<typename Src, typename Boxes>
void ParticlesExporter<AoSCMTS, CPU>::delete_particles_in(Src& src, Boxes const& boxes)
{
    src.erase_if([&](auto const& p) { return isIn(p.iCell(), boxes); });
}

template<>
template<typename Dst, typename Src, std::size_t dim>
void ParticlesExporter<AoSCMTS, CPU>::move_in_domain(Dst& dst, Src& src,
                                                     Box<int, dim> const& domain_box)
{
    src.erase_if([&](auto const& p) { return isIn(p.iCell(), domain_box); },
                 [&](auto const& p) { dst.push_back(p); });
    dst.template on_appended<ParticleType::Domain>();
}

template<>
template<typename Dst, typename Src, std::size_t dim>
void ParticlesExporter<AoSCMTS, CPU>::move_in_ghost_layer(Dst& dst, Src& src,
                                                          Box<int, dim> const& domain_box,
                                                          Box<int, dim> const& ghost_box)
{
    auto const ghost_layer_boxes = ghost_box.remove(domain_box);
    src.erase_if([&](auto const& p) { return isIn(p.iCell(), ghost_layer_boxes); },
                 [&](auto const& p) { dst.push_back(p); });
    dst.template on_appended<ParticleType::Domain>();
}

template<>
template<typename Dst, typename Src, std::size_t dim, typename Boxes>
void ParticlesExporter<AoSCMTS, CPU>::move_in_ghost_layer(Dst& dst, Src& src,
                                                          Box<int, dim> const& domain_box,
                                                          Boxes const& ghost_boxes)
{
    for (auto const& gb : ghost_boxes)
        this->move_in_ghost_layer(dst, src, domain_box, gb);
}


} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_EXPORTER */
