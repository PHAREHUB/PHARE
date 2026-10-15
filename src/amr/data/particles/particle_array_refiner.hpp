// to #include "core/data/particles/particle_array_exporter.hpp"

#ifndef PHARE_AMR_DATA_PARTICLES_PARTICLE_ARRAY_REFINER
#define PHARE_AMR_DATA_PARTICLES_PARTICLE_ARRAY_REFINER

#include "core/data/particles/particle_array_def.hpp"
#include "amr/data/particles/refining/particles_refining.hpp"

namespace PHARE::amr
{


template<auto options, auto type = core::ParticleType::Domain, typename Src, typename Dst>
void export_refined_particles(Src const& src, Dst& dst, auto const& dst_boxes, auto refiner,
                              auto transformer, auto const& dst_amr_box)
{
    using Refiner = ParticlesRefiner<Src::layout_mode, Src::alloc_mode>;

    std::string_view constexpr static FN_ID     = "export_particles,";
    [[maybe_unused]] auto constexpr function_id = core::join_string_views_v<FN_ID, Src::type_id>;
    PHARE_LOG_SCOPE(3, function_id);

    RefineArgs<Src, Dst> args{dst, src, dst_boxes, dst_amr_box};
    Refiner{}.template operator()<options, type>(args, refiner, transformer);
}

} // namespace PHARE::amr

#endif /* PHARE_AMR_DATA_PARTICLES_PARTICLE_ARRAY_REFINER */
