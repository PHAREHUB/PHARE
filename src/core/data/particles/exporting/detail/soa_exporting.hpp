#ifndef PHARE_CORE_DATA_PARTICLES_EXPORTING_DETAIL_SOA_EXPORTER
#define PHARE_CORE_DATA_PARTICLES_EXPORTING_DETAIL_SOA_EXPORTER


#include "core/utilities/box/box.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/exporting/detail/def_exporting.hpp"

#include <cstddef>


namespace PHARE::core::detail::soa_exporting
{

template<typename Src>
void erase_if(Src& src, auto const& erase, auto const& on_erase)
{
    std::size_t n = src.size();
    for (std::size_t i = 0; i < n;)
    {
        if (erase(src.iCell(i)))
        {
            on_erase(i);
            src.swap(i, --n);
        }
        else
            ++i;
    }
    src.resize(n);
}

template<typename Src>
void erase_if(Src& src, auto const& erase)
{
    erase_if(src, erase, [](auto) {});
}

template<typename Dst, typename Src>
void move_if(Dst& dst, Src& src, auto const& move)
{
    erase_if(src, move, [&](auto const i) {
        if constexpr (Dst::layout_mode == LayoutMode::SoAVX)
            dst.push_back(src.copy(i));
        else
            dst.emplace_back(src.copy(i));
    });
}

} // namespace PHARE::core::detail::soa_exporting


namespace PHARE::core
{

using enum LayoutMode;
using enum AllocatorMode;


#define PHARE_SOA_EXPORTER(LAYOUT)                                                                 \
    template<>                                                                                     \
    template<typename Src, std::size_t dim>                                                        \
    void ParticlesExporter<LAYOUT, CPU>::delete_particles_in(Src& src, Box<int, dim> const& box)   \
    {                                                                                              \
        detail::soa_exporting::erase_if(src, [&](auto const& cell) { return isIn(cell, box); });   \
    }                                                                                              \
                                                                                                   \
    template<>                                                                                     \
    template<typename Src, typename Boxes>                                                         \
    void ParticlesExporter<LAYOUT, CPU>::delete_particles_in(Src& src, Boxes const& boxes)         \
    {                                                                                              \
        detail::soa_exporting::erase_if(src, [&](auto const& cell) { return isIn(cell, boxes); }); \
    }                                                                                              \
                                                                                                   \
    template<>                                                                                     \
    template<typename Src, std::size_t dim>                                                        \
    void ParticlesExporter<LAYOUT, CPU>::delete_particles_not_in(Src& src,                         \
                                                                 Box<int, dim> const& box)         \
    {                                                                                              \
        detail::soa_exporting::erase_if(src, [&](auto const& cell) { return !isIn(cell, box); });  \
    }                                                                                              \
                                                                                                   \
    template<>                                                                                     \
    template<typename Src, typename Boxes>                                                         \
    void ParticlesExporter<LAYOUT, CPU>::delete_particles_not_in(Src& src, Boxes const& boxes)     \
    {                                                                                              \
        detail::soa_exporting::erase_if(src,                                                       \
                                        [&](auto const& cell) { return !isIn(cell, boxes); });     \
    }                                                                                              \
                                                                                                   \
    template<>                                                                                     \
    template<typename Dst, typename Src, std::size_t dim>                                          \
    void ParticlesExporter<LAYOUT, CPU>::move_in_domain(Dst& dst, Src& src,                        \
                                                        Box<int, dim> const& domain_box)           \
    {                                                                                              \
        detail::soa_exporting::move_if(dst, src,                                                   \
                                       [&](auto const& cell) { return isIn(cell, domain_box); });  \
    }                                                                                              \
                                                                                                   \
    template<>                                                                                     \
    template<typename Dst, typename Src, std::size_t dim>                                          \
    void ParticlesExporter<LAYOUT, CPU>::move_in_ghost_layer(                                      \
        Dst& dst, Src& src, Box<int, dim> const& domain_box, Box<int, dim> const& ghost_box)       \
    {                                                                                              \
        auto const ghost_layer_boxes = ghost_box.remove(domain_box);                               \
        detail::soa_exporting::move_if(                                                            \
            dst, src, [&](auto const& cell) { return isIn(cell, ghost_layer_boxes); });            \
    }                                                                                              \
                                                                                                   \
    template<>                                                                                     \
    template<typename Dst, typename Src, std::size_t dim, typename Boxes>                          \
    void ParticlesExporter<LAYOUT, CPU>::move_in_ghost_layer(                                      \
        Dst& dst, Src& src, Box<int, dim> const& domain_box, Boxes const& ghost_boxes)             \
    {                                                                                              \
        for (auto const& gb : ghost_boxes)                                                         \
            this->move_in_ghost_layer(dst, src, domain_box, gb);                                   \
    }

PHARE_SOA_EXPORTER(SoA)
PHARE_SOA_EXPORTER(SoAVX)

#undef PHARE_SOA_EXPORTER

} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_EXPORTING_DETAIL_SOA_EXPORTER */
