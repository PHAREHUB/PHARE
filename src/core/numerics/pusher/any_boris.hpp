#ifndef PHARE_CORE_PUSHER_BORIS_SIMPLER_HPP
#define PHARE_CORE_PUSHER_BORIS_SIMPLER_HPP

#include "core/utilities/kernels.hpp"
#include "core/numerics/pusher/boris/basics.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_array_service.hpp"


#include <cmath>
#include <array>
#include <cassert>
#include <cstddef>
#include <iterator>
#include <algorithm>
#include <stdexcept>


namespace PHARE::core
{

template<LayoutMode layout, AllocatorMode alloc, std::size_t dim,
         typename Interpolator, typename GridLayout>
struct AnyBorisImpl; // defined in boris/mkn/mkn_any_boris.hpp (included at bottom)


template<std::size_t dim, typename Interpolator, typename GridLayout>
class AnyBorisPusher
{
    auto static mesh(std::array<double, dim> const& ms, double const& ts)
    {
        std::array<double, dim> halfDtOverDl;
        std::transform(std::begin(ms), std::end(ms), std::begin(halfDtOverDl),
                       [ts](auto const& x) { return 0.5 * ts / x; });
        return halfDtOverDl;
    }

public:
    AnyBorisPusher(                        //
        GridLayout const& layout,          //
        std::array<double, dim> const& ms, //
        double const& ts,                  //
        double const& mass)
        : layout_{layout}
        , halfDtOverDl_{mesh(ms, ts)}
        , mass_{mass}
        , dt_{ts}
    {
        PHARE_ASSERT(mass_ > 0);
        PHARE_ASSERT(dt_ > 0);
        PHARE_ASSERT(dto2m > 0);
    }


    template<typename Particles, typename Electromag>
    void operator()(Particles& particles, Electromag const& em) const
    {
        move<ParticleType::Domain>(particles, em);
    }

    template<auto type, typename Particles, typename Electromag>
    void push(Particles& particles, Electromag const& em) const
    {
        move<type>(particles, em);
    }


    template<auto type, typename Particles, typename Electromag>
    void move(Particles& particles, Electromag const& em) const _PHARE_ALL_FN_
    {
        if (particles.size() == 0)
            return;

        PHARE_LOG_SCOPE(1, "boris::move");

        using Impl = AnyBorisImpl<Particles::layout_mode, Particles::alloc_mode,
                                   dim, Interpolator, GridLayout>;
        Impl::template move<type>(particles, em, interpolator_, layout_, halfDtOverDl_, dto2m);
    }


    mutable Interpolator interpolator_;
    GridLayout const& layout_;
    std::array<double, dim> const halfDtOverDl_;
    double const mass_ = 0;
    double const dt_   = 0;
    double const dto2m = 0.5 * dt_ / mass_;
};

} // namespace PHARE::core


#include "core/numerics/pusher/boris/mkn/mkn_any_boris.hpp"


#endif /* PHARE_CORE_PUSHER_BORIS_SIMPLER_HPP */
