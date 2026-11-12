#ifndef PHARE_AMR_DATA_PARTICLES_PARTICLE_ARRAY_REFINER
#define PHARE_AMR_DATA_PARTICLES_PARTICLE_ARRAY_REFINER

#include "core/data/particles/particle_array.hpp"
#include "core/utilities/box/box_spans.hpp"
#include "core/data/particles/particle_array_selector.hpp"

#include "amr/utilities/box/amr_box.hpp"

#include <SAMRAI/hier/BoxContainer.h>

#include <vector>

namespace PHARE::amr
{

// --- public API ---

// dst_amr_box: destination patch's interior box, distinct from dst_boxes - lets the Ghost
// reserve estimate exclude particles whose offspring only ever land in the interior.
template<typename Src, typename Dst>
struct RefineArgs
{
    using Box_t = core::Box<int, Src::dimension>;

    Dst& dst;
    Src const& src;
    SAMRAI::hier::BoxContainer const& dst_boxes;
    Box_t const& dst_amr_box;
};

// hides the per layout/alloc impl behind operator() - see implementation section below
template<auto layout_mde, auto alloc_mde>
struct ParticlesRefiner;


template<auto options, auto type = core::ParticleType::Domain, typename Src, typename Dst>
void export_refined_particles(Src const& src, Dst& dst, auto const& dst_boxes, auto refiner,
                              auto transformer, auto const& dst_amr_box)
{
    using Refiner = ParticlesRefiner<Src::layout_mode, Src::alloc_mode>;

    RefineArgs<Src, Dst> args{dst, src, dst_boxes, dst_amr_box};
    Refiner{}.template operator()<options, type>(args, refiner, transformer);
}


// --- implementations ---

template<auto layout_mde, auto alloc_mde>
struct ParticlesRefiner
{
    static_assert(core::all_are<core::LayoutMode>(layout_mde));
    static_assert(core::all_are<core::AllocatorMode>(alloc_mde));

    auto constexpr static layout_mode = layout_mde;
    auto constexpr static alloc_mode  = alloc_mde;

    template<auto options, auto type, typename Src, typename Dst>
    void operator()(RefineArgs<Src, Dst>&, auto /*Refiner*/, auto /*Transformer*/);
};


using enum core::LayoutMode;
using enum core::AllocatorMode;


// for Ghost, particles whose offspring only ever land in the interior dst_amr_box never
// reach the ghost shell, so exclude them via inclusion-exclusion.
template<auto type, typename Src, typename Box_t>
std::size_t relevant_coarse_count(Src const& src, Box_t const& coarseDstBox,
                                  Box_t const& dst_amr_box)
{
    auto const total = core::count_particles(src, coarseDstBox);
    if constexpr (type == core::ParticleType::Ghost)
    {
        auto const coarseAmrBox = coarsen_box(dst_amr_box);
        if (auto const overlap = coarseDstBox * coarseAmrBox)
            return total - core::count_particles(src, *overlap);
    }
    return total;
}

// every one of count's nbRefinedPart children lands somewhere in the box coarseDstBox was
// derived from (that's the whole box, not one fine cell) - no /2**dim here, that would turn
// a box total into a per-cell average.
template<auto options, typename Src>
constexpr std::size_t scale_by_density(std::size_t const count, double const factor)
{
    return static_cast<std::size_t>(count * options.nbRefinedPart * factor);
}

// per-cell reserve factor for a coarse cell's 2**dim fine children: the even spread average
template<auto options, auto type, typename Src>
constexpr double per_cell_density_factor()
{
    return 1. / (std::size_t{1} << Src::dimension);
}

template<auto type, auto options, typename Src, typename Box_t>
std::size_t reserve_flat_count(Src const& src, Box_t const& dst_box, Box_t const& dst_amr_box,
                               double const factor)
{
    // Domain children are clipped to the exact domain box during streaming
    // (ParticlesRefining::forBoxes' per_particle partitions against domainBox). Level-ghost
    // overlaps always sit on top of coarser *domain* cells (AMR nesting guarantees a fine
    // level's ghost region is covered by the coarse level's interior, never its ghost layer),
    // so neither case needs the particle_ghost_width margin here.
    auto const coarseDstBox = coarsen_box(dst_box);
    auto const count        = relevant_coarse_count<type>(src, coarseDstBox, dst_amr_box);
    return scale_by_density<options, Src>(count, factor);
}

template<auto type, auto options, typename Src, typename Dst, typename Box_t>
void reserve_flat(Src const& src, Dst& dst, Box_t const& dst_box, Box_t const& dst_amr_box,
                  double const factor)
{
    dst.reserve(dst.size() + reserve_flat_count<type, options>(src, dst_box, dst_amr_box, factor));
}

template<auto type, auto options, typename Src, typename Dst, typename Box_t>
void reserve_flat(Src const& src, Dst& dst, SAMRAI::hier::BoxContainer const& dst_boxes,
                  Box_t const& dst_amr_box)
{
    auto constexpr domain_fac_per_layout = std::array{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    std::size_t n = 0;
    for (auto const& samrai_box : dst_boxes)
    {
        n += reserve_flat_count<type, options>(src, phare_box_from<Src::dimension>(samrai_box),
                                               dst_amr_box, .66);

        // grow box by particle_ghost_width
        // remove samrai_box from new bigger box
        // operate across remainder at lower weight
        if constexpr (type != core::ParticleType::Domain)
        {
            SAMRAI::hier::BoxContainer ghost_boxes{grow(samrai_box, options.particle_ghost_width)};
            ghost_boxes.removeIntersections(samrai_box);
            for (auto const& border_box : ghost_boxes)
                n += reserve_flat_count<type, options>(
                    src, phare_box_from<Src::dimension>(border_box), dst_amr_box, .1);
        }
    }

    dst.reserve(dst.size() + n);
}

template<auto type, auto options, typename Src, typename Dst, typename Box_t>
void stream_split_particles(Src const& src, Dst& dst, Box_t const& dst_box, auto fn0, auto fn1)
{
    std::uint16_t constexpr static N       = 256;
    static constexpr auto base_layout_type = core::base_layout_type<Src>();
    static constexpr auto array_opts
        = Src::options.with_storage(core::StorageMode::ARRAY).with_layout(base_layout_type);
    static constexpr auto array_type_opts
        = core::ParticleArrayTypeOptions<array_opts, base_layout_type,
                                         core::StorageMode::ARRAY>::FROM(Src::options, N);
    using ArrayParticleArray = core::ParticleArrayResolver<array_opts, array_type_opts>::value_type;

    // Resolved the same way as ArrayParticleArray above (not Src::Span_t): that alias
    // only exists on the flat AoS/SoA/SoAVX resolved classes.
    static constexpr auto span_opts
        = Src::options.with_storage(core::StorageMode::SPAN).with_layout(base_layout_type);
    static constexpr auto span_type_opts
        = core::ParticleArrayTypeOptions<span_opts, base_layout_type,
                                         core::StorageMode::SPAN>::FROM(Src::options);
    using SpanParticleArray = core::ParticleArrayResolver<span_opts, span_type_opts>::value_type;

    // particle_ghost_width == Splitter::maxCellDistanceFromSplit() (see splitter.hpp) - a
    // coarse particle up to that many fine cells outside dst_box can still have a child land
    // inside it, so this margin is load-bearing for correctness, unlike reserve_flat_count's.
    auto const splitBox     = grow(dst_box, options.particle_ghost_width);
    auto const coarseDstBox = coarsen_box(splitBox);

    std::uint16_t big_buffer_cnt = 0;
    ArrayParticleArray big_buffer;
    SpanParticleArray span{big_buffer};

    auto const send = [&]() {
        assert(big_buffer_cnt <= N);
        span.resize(big_buffer_cnt);
        append_particles<type>(span, dst);
        big_buffer_cnt = 0;
    };

    core::per_particle(src, [&](auto const& particle) {
        if (not isIn(particle, coarseDstBox))
            return;

        if (not isIn(fn0(particle), splitBox))
            return;

        auto&& [p_count, buffer] = fn1(particle, dst_box);
        if (big_buffer_cnt + p_count > N)
            send();

        std::copy(buffer.data(), buffer.data() + p_count, big_buffer.data() + big_buffer_cnt);
        big_buffer_cnt += p_count;
    });

    if (big_buffer_cnt)
        send();
}


template<>
template<auto options, auto type, typename Src, typename Dst>
void ParticlesRefiner<AoS, CPU>::operator()(RefineArgs<Src, Dst>& args, auto fn0, auto fn1)
{
    auto&& [dst, src, dst_boxes, dst_amr_box] = args;
    reserve_flat<type, options>(src, dst, dst_boxes, dst_amr_box);
    for (auto const& samrai_box : dst_boxes)
    {
        auto const dst_box = phare_box_from<Src::dimension>(samrai_box);
        stream_split_particles<type, options>(src, dst, dst_box, fn0, fn1);
    }
}

template<>
template<auto options, auto type, typename Src, typename Dst>
void ParticlesRefiner<AoSMapped, CPU>::operator()(RefineArgs<Src, Dst>& args, auto fn0, auto fn1)
{
    ParticlesRefiner<AoS, CPU>{}.template operator()<options, type>(args, fn0, fn1);
}


template<auto type, typename Dst, typename Tile, typename Box_t>
std::optional<Box_t> refine_tile_overlap(Dst const& dst, Tile const& dst_tile, Box_t const& dst_box,
                                         std::size_t const ghost_nbr)
{
    // dst_box can extend past this array's own domain (periodic wrap, a restricted/
    // temporary ParticlesData SAMRAI builds internally, or the level ghost layer).
    // Particles are never duplicated across tiles, so each bound's margin is owned
    // by whichever tile touches the domain edge it's adjacent to (same clamp rule
    // tag_cells_ uses, so TileSet::at(ghost_cell) finds them); a tile that doesn't
    // own a bound just clips to its own edge there instead of the whole region
    // being rejected. Readers needing a neighbour's level ghosts gather via at().
    auto const growbox = grow(dst_tile, ghost_nbr);
    if (not(dst_box * growbox))
        return std::nullopt;

    auto const& domain_box = dst.box();
    auto const& tile_box   = *dst_tile;
    auto region            = dst_box;
    for (std::size_t d = 0; d < Dst::dimension; ++d)
    {
        if (dst_box.lower[d] < domain_box.lower[d])
        {
            bool const owns_edge = domain_box.lower[d] >= tile_box.lower[d]
                                   and domain_box.lower[d] <= tile_box.upper[d];
            region.lower[d]
                = owns_edge ? dst_box.lower[d] : std::max(dst_box.lower[d], tile_box.lower[d]);
        }
        else
            region.lower[d] = std::max(dst_box.lower[d], tile_box.lower[d]);

        if (dst_box.upper[d] > domain_box.upper[d])
        {
            bool const owns_edge = domain_box.upper[d] >= tile_box.lower[d]
                                   and domain_box.upper[d] <= tile_box.upper[d];
            region.upper[d]
                = owns_edge ? dst_box.upper[d] : std::min(dst_box.upper[d], tile_box.upper[d]);
        }
        else
            region.upper[d] = std::min(dst_box.upper[d], tile_box.upper[d]);

        if (region.lower[d] > region.upper[d])
            return std::nullopt;
    }
    return region;
}



template<>
template<auto options, auto type, typename Src, typename Dst>
void ParticlesRefiner<AoSCMTS, CPU>::operator()(RefineArgs<Src, Dst>& args, auto fn0, auto fn1)
{
    // one flat cell-mapped array per tile: each dst tile has a single flat buffer, so
    // reserve per tile overlap
    auto&& [dst, src, dst_boxes, dst_amr_box] = args;
    auto const ghost_nbr                      = options.particle_ghost_width;

    static_assert(type == core::ParticleType::Domain || type == core::ParticleType::Ghost);

    using Box_t     = core::Box<int, Src::dimension>;
    using DstTile_t = std::decay_t<decltype(dst()[0])>;
    struct DstOverlap
    {
        DstTile_t* tile;
        Box_t overlap;      // fine cells of dst_box this dst tile receives
        Box_t coarse_reach; // coarse cells whose particles can split into overlap
    };
    thread_local std::vector<DstOverlap> dst_overlaps; // reused, no per call allocation

    for (auto const& samrai_box : dst_boxes)
    {
        auto const dst_box      = phare_box_from<Src::dimension>(samrai_box);
        auto const coarseDstBox = coarsen_box(dst_box);

        dst_overlaps.clear();
        for (auto& dst_tile : dst())
            if (auto const overlap = refine_tile_overlap<type>(dst, dst_tile, dst_box, ghost_nbr))
            {
                reserve_flat<type, options>(src, dst_tile(), *overlap, dst_amr_box,
                                            per_cell_density_factor<options, type, Src>());
                // same coarse box stream_split_particles filters source particles with
                dst_overlaps.emplace_back(&dst_tile, *overlap,
                                          coarsen_box(grow(*overlap, ghost_nbr)));
            }

        for (auto const& src_tile : src())
        {
            auto const src_reach = grow(src_tile, ghost_nbr); // cells src_tile() may hold
            if (not(coarseDstBox * src_reach))
                continue;

            // only dst tiles whose coarse reach touches src_tile can receive its children,
            // skip scanning src_tile's particles for every other dst tile
            for (auto const& [dst_tile, overlap, coarse_reach] : dst_overlaps)
                if (coarse_reach * src_reach)
                    stream_split_particles<type, options>(src_tile(), (*dst_tile)(), overlap,
                                                          fn0, fn1);
        }
    }

    dst.template on_appended<type>();
}


} // namespace PHARE::amr

#endif /* PHARE_AMR_DATA_PARTICLES_PARTICLE_ARRAY_REFINER */
