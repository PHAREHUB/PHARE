// to #include "core/data/particles/particle_array_exporter.hpp"

#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_EXPORTER
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_EXPORTER

#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/exporting/particles_exporting.hpp"
#include "core/data/particles/exporting/detail/def_exporting.hpp"

namespace PHARE::core
{

template<typename Src, typename Box_t>
void delete_particles_in(Src& src, Box_t const& box)
{
    using Exporter = ParticlesExporter<Src::layout_mode, Src::alloc_mode>;

    std::string_view constexpr static FN_ID     = "delete_particles,";
    [[maybe_unused]] auto constexpr function_id = join_string_views_v<FN_ID, Src::type_id>;
    PHARE_LOG_SCOPE(3, function_id);

    Exporter{}.delete_particles_in(src, box);
}

template<typename Src, typename Box_t>
void delete_particles_not_in(Src& src, Box_t const& box)
{
    using Exporter = ParticlesExporter<Src::layout_mode, Src::alloc_mode>;

    std::string_view constexpr static FN_ID     = "delete_particles,";
    [[maybe_unused]] auto constexpr function_id = join_string_views_v<FN_ID, Src::type_id>;
    PHARE_LOG_SCOPE(3, function_id);

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


} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_EXPORTER */
