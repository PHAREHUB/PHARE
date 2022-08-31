#ifndef PHARE_CORE_NUMERICS_PUSHER_PUSHER_HPP
#define PHARE_CORE_NUMERICS_PUSHER_PUSHER_HPP

#include <cstddef>
#include <type_traits>
#include <utility>
#include <functional>

#include "core/utilities/range/range.hpp"
#include "core/data/particles/particle.hpp"

namespace PHARE
{
namespace core
{
    template<std::size_t dim, typename ParticleRange, typename Electromag, typename Interpolator,
             typename BoundaryCondition, typename GridLayout, bool UseParticleIndex = true>
    class Pusher
    {
    protected:
        static auto constexpr dimension = GridLayout::dimension;

        using ParticleSelector = std::function<ParticleRange(ParticleRange&)>;

    public:
        /**
         * @brief setMeshAndTimeStep allows to let the pusher know what is the mesh
         * size and time step in the domain where particles are to be pushed.
         */
        virtual void setMeshAndTimeStep(std::array<double, dim> const& ms,
                                        double const ts) _PHARE_ALL_FN_ = 0;


        // TODO : to really be independant on boris which has 2 push steps
        // we should have an arbitrary number of selectors, 1 per push step
        virtual ParticleRange move(ParticleRange const& rangeIn, ParticleRange& rangeOut,
                                   Electromag const& emFields, double mass,
                                   Interpolator& interpolator, GridLayout const& layout,
                                   ParticleSelector firstSelector, ParticleSelector secondSelector)
            = 0;



        Pusher() _PHARE_ALL_FN_ {}
        virtual ~Pusher() _PHARE_ALL_FN_ {}
    };

} // namespace core
} // namespace PHARE

#endif
