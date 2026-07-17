// to #include "core/data/particles/particle_array_exporter.hpp"

#ifndef PHARE_AMR_DATA_PARTICLES_PARTICLE_ARRAY_REFINER
#define PHARE_AMR_DATA_PARTICLES_PARTICLE_ARRAY_REFINER

#include "amr/data/particles/refining/detail/def_refining.hpp"
#include "amr/data/particles/refining/particles_refining.hpp"

namespace PHARE::amr
{


template<auto type, typename Src, typename Dst, typename Box_t, typename RefineFn,
         typename Transformer>
void export_refined_particles( //
    Src const& src, Dst& dst, Box_t const& box, RefineFn refiner, Transformer fn0)
{
    using Refiner = ParticlesRefiner<Src::layout_mode, Src::alloc_mode>;

    std::string_view constexpr static FN_ID     = "export_particles,";
    [[maybe_unused]] auto constexpr function_id = core::join_string_views_v<FN_ID, Src::type_id>;
    PHARE_LOG_SCOPE(3, function_id);

    Refiner{}.template operator()<type>(src, dst, box, refiner, fn0);
}

} // namespace PHARE::amr

#endif /* PHARE_AMR_DATA_PARTICLES_PARTICLE_ARRAY_REFINER */
