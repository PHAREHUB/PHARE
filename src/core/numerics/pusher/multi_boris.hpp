#ifndef PHARE_CORE_PUSHER_MULTI_BORIS_HPP
#define PHARE_CORE_PUSHER_MULTI_BORIS_HPP

#if PHARE_HAVE_MKN_GPU

#include "core/def.hpp"
#include "core/utilities/span.hpp"
#include "core/utilities/kernels.hpp"
#include "core/utilities/thread_pool.hpp"
#include "core/data/electromag/electromag.hpp"
#include "core/data/tensorfield/tensorfield.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_array_service.hpp"
#include "core/numerics/pusher/boris/basics.hpp"
#include "core/numerics/interpolator/interpolating.hpp"

// MultiBoris state struct + MultiBorisFunctors + MultiBorisPusherImpl<AoSTS/AoSMapped/AoSPC>
#include "core/numerics/pusher/boris/mkn/mkn_multi_boris.hpp"

#include <cmath>
#include <memory>
#include <cassert>
#include <cstddef>
#include <functional>


namespace PHARE::core
{

// Thin dispatcher: routes to MultiBorisPusherImpl<layout_mode, alloc_mode, ...>.
template<typename GridLayout, typename Particles, typename Electromag, typename Interpolator>
class MultiBorisPusher
{
    static constexpr auto layout = Particles::layout_mode;
    static constexpr auto alloc  = Particles::alloc_mode;
    using Impl
        = MultiBorisPusherImpl<layout, alloc, GridLayout, Particles, Electromag, Interpolator>;

public:
    template<MultiBorisMode mode = MultiBorisMode::REF, typename ModelAccessor>
    static void move(MultiBoris<ModelAccessor>& in, auto const& boxings)
    {
        Impl::template move<mode>(in, boxings);
    }
};

} // namespace PHARE::core

#endif // PHARE_HAVE_MKN_GPU

#endif /* PHARE_CORE_PUSHER_MULTI_BORIS_HPP */
