#ifndef PHARE_CORE_NUMERICS_PUSHER_PUSHER_H
#define PHARE_CORE_NUMERICS_PUSHER_PUSHER_H

#include <cstddef>
#include <utility>
#include <functional>

#include "core/utilities/range/range.h"
#include "core/data/particles/particle.h"

namespace PHARE
{
namespace core
{
    template<std::size_t dim, typename ParticleIterator, typename Electromag, typename Interpolator,
             typename BoundaryCondition, typename GridLayout>
    class Pusher
    {
    protected:
        using ParticleRange             = Range<ParticleIterator>;
        static auto constexpr dimension = GridLayout::dimension;
        using ParticleSelector          = std::function<bool(Particle<dimension> const&)>;

    public:
        /** Move all particles in rangeIn from t=n to t=n+1 and store their new
         * position in rangeOut.
         *
         * Particles in rangeOut are placed according to the return of the ParticleSelector.
         * Particles for which the selector retuns true are placed before those for which it
         * is false. The pivot is the iterator separating the two parts. Here it is assumed
         * the selector returns true for particles staying in the patch and false otherwise.
         *
         * move() assumes that particles in rangeIn and those in rangeOut have the same
         * state upon entering the function.
         *
         * @param rangeIn : range of iterators on particles at time t=n to be pushed
         * @param rangeOut: output range of iterators on particles at t=n+1. rangeOut
         * must have the same size as rangeIn (nbrParticles(rangeIn) == nbrParticles(rangeOut).
         * @param E: electric vector field used to accelerate particles
         * @param B: magnetic vector field used to accelerate particles
         * @param selector : used to place particles in rangeOut.
         * @param bc : physical boundary condition. Manage particles that intersect with a
         * physical domain bounday.
         *
         * @return the function returns an ParticleArray::iterator on the first particle
         * for which the selector returns false. The selector returns true for Particles
         * in range [rangeOut.begin, pivot[, and false for those in [pivot,rangeOut.end[
         */
        virtual ParticleIterator
        move(ParticleRange const& rangeIn, ParticleRange& rangeOut, Electromag const& emFields,
             double mass, Interpolator& interpolator, ParticleSelector const& particleIsNotLeaving,
             BoundaryCondition& bc, GridLayout const& layout)
            = 0;


        /** this overload of move() is used for particles for which one knows that
         * they will not need a boundary condition treatment. Particles in rangeOut
         * are sorted according to whether they are detected as leaving
         *
         * This overload is typically used to push particles outside the domain, like
         * ghost particles.
         */
        virtual ParticleIterator // decltype(std::declval<ParticleRange>().end())
        move(ParticleRange const& rangeIn, ParticleRange& rangeOut, Electromag const& emFields,
             double mass, Interpolator& interpolator, ParticleSelector const& particleIsNotLeaving,
             GridLayout const& layout)
            = 0;


        virtual ParticleIterator move(ParticleRange& range, Electromag const& emFields, double mass,
                                      Interpolator& interpolator,
                                      ParticleSelector const& particleIsNotLeaving,
                                      GridLayout const& layout)
            = 0;


        /**
         * @brief setMeshAndTimeStep allows to let the pusher know what is the mesh
         * size and time step in the domain where particles are to be pushed.
         */
        virtual void setMeshAndTimeStep(std::array<double, dim> const& ms,
                                        double ts) _PHARE_FN_SIG_ = 0;

        Pusher() _PHARE_FN_SIG_ {}
        virtual ~Pusher() _PHARE_FN_SIG_ {}
    };

} // namespace core
} // namespace PHARE

#endif
