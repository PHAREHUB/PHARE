#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_APPENDER
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_APPENDER

#include "core/data/particles/appending/particles_appending.hpp"

namespace PHARE::core
{


template<auto type, typename Src, typename Dst>
void append_particles(Src const& src, Dst& dst)
{
    using Appending
        = ParticlesAppender<Src::layout_mode, Src::alloc_mode, Dst::layout_mode, Dst::alloc_mode>;

    PHARE_DEBUG_DO(int const old_size = dst.size();)

    std::string_view constexpr static FN_ID     = "append_particles,";
    [[maybe_unused]] auto constexpr function_id = join_string_views_v<FN_ID, Dst::type_id>;
    PHARE_LOG_SCOPE(3, function_id);

    Appending{0, src.size()}.template operator()<type>(src, dst);

    assert(dst.size() == old_size + src.size());
}


} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_APPENDER */
