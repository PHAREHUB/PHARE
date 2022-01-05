#ifndef PHARE_CORE_PUSHER_BORIS_H
#define PHARE_CORE_PUSHER_BORIS_H

#include <array>
#include <cmath>
#include <cstddef>
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include "core/numerics/pusher/pusher.hpp"
#include "core/utilities/range/range.hpp"
#include "core/errors.hpp"
#include "core/logger.hpp"
#include "core/data/particles/particle.hpp"

namespace PHARE
{
namespace core
{
    template<std::size_t dim, typename ParticleIterator, typename Particle_t, typename Electromag,
             typename Interpolator, typename BoundaryCondition, typename GridLayout>
    class BorisPusher : public Pusher<dim, ParticleIterator, Particle_t, Electromag, Interpolator,
                                      BoundaryCondition, GridLayout>
    {
    public:
        using Super            = Pusher<dim, ParticleIterator, Particle_t, Electromag, Interpolator,
                             BoundaryCondition, GridLayout>;
        using ParticleSelector = typename Super::ParticleSelector;
        using ParticleRange    = Range<ParticleIterator>;

    private:
    public:
        // This move function should be considered when being used so that all particles are pushed
        // twice - see: https://github.com/PHAREHUB/PHARE/issues/571
        /** see Pusher::move() documentation*/
        ParticleIterator move(ParticleRange const& rangeIn, ParticleRange& rangeOut,
                              Electromag const& emFields, double mass, Interpolator& interpolator,
                              ParticleSelector const& particleIsNotLeaving, BoundaryCondition& bc,
                              GridLayout const& layout) override
        {
            // push the particles of half a step
            // rangeIn : t=n, rangeOut : t=n+1/Z
            // get a pointer on the first particle of rangeOut that leaves the patch
            auto firstLeaving
                = pushStep_(rangeIn, rangeOut, particleIsNotLeaving, PushStep::PrePush);

            // apply boundary condition on the particles in [firstLeaving, rangeOut.end[
            // that actually leave through a physical boundary condition
            // get a pointer on the new end of rangeOut. Particles passed newEnd
            // are those that have left the patch through a non-physical boundary
            // they should be discarded now
            auto newEnd = bc.applyOutgoingParticleBC(firstLeaving, rangeOut.end());

            rangeOut = makeRange(rangeOut.begin(), std::move(newEnd));

            // get electromagnetic fields interpolated on the particles of rangeOut
            // stop at newEnd.
            interpolator(rangeOut.begin(), rangeOut.end(), emFields, layout);

            // get the particle velocity from t=n to t=n+1
            accelerate_(rangeOut, rangeOut, mass);

            // now advance the particles from t=n+1/2 to t=n+1 using v_{n+1} just calculated
            // and get a pointer to the first leaving particle
            firstLeaving = pushStep_(rangeOut, rangeOut, particleIsNotLeaving, PushStep::PostPush);

            // apply BC on the leaving particles that leave through physical BC
            // and get pointer on new End, discarding particles leaving elsewhere
            newEnd = bc.applyOutgoingParticleBC(firstLeaving, rangeOut.end());

            rangeOut = makeRange(rangeOut.begin(), std::move(newEnd));

            return rangeOut.end();
        }


        /** see Pusher::move() documentation*/
        ParticleIterator move(ParticleRange const& rangeIn, ParticleRange& rangeOut,
                              Electromag const& emFields, double mass, Interpolator& interpolator,
                              ParticleSelector const& particleIsNotLeaving,
                              GridLayout const& layout) override
        {
            PHARE_LOG_SCOPE("Boris::move_no_bc");

            // push the particles of half a step
            // rangeIn : t=n, rangeOut : t=n+1/2
            // Do not partition on this step - this is to keep all domain and ghost
            //   particles consistent. see: https://github.com/PHAREHUB/PHARE/issues/571
            pushStep_(rangeIn, rangeOut, PushStep::PrePush);

            //  get electromagnetic fields interpolated on the particles of rangeOut
            //  stop at newEnd.
            interpolator(rangeOut.begin(), rangeOut.end(), emFields, layout);

            //  get the particle velocity from t=n to t=n+1
            accelerate_(rangeOut, rangeOut, mass);

            // now advance the particles from t=n+1/2 to t=n+1 using v_{n+1} just calculated
            // and get a pointer to the first leaving particle
            auto firstLeaving
                = pushStep_(rangeOut, rangeOut, particleIsNotLeaving, PushStep::PostPush);

            rangeOut = makeRange(rangeOut.begin(), std::move(firstLeaving));

            return rangeOut.end();
        }

        /** see Pusher::move() documentation*/
        // This function used used during particle ghost push, so that any particles outside the
        //  ghostbox after first pushStep are ignored, and only particles in the ghostLayer are
        //  kept after second pushStep
        ParticleIterator move(ParticleRange const& rangeIn, ParticleRange& rangeOut,
                              Electromag const& emFields, double mass, Interpolator& interpolator,
                              ParticleSelector const& particleIsNotLeaving0,
                              ParticleSelector const& particleIsNotLeaving1,
                              GridLayout const& layout) override
        {
            PHARE_LOG_SCOPE("Boris::move_no_bc");

            // push the particles of half a step
            // rangeIn : t=n, rangeOut : t=n+1/2
            auto firstLeaving
                = pushStep_(rangeIn, rangeOut, particleIsNotLeaving0, PushStep::PrePush);

            rangeOut = makeRange(rangeOut.begin(), std::move(firstLeaving));

            // get electromagnetic fields interpolated on the particles of rangeOut stop at newEnd.
            interpolator(rangeOut.begin(), rangeOut.end(), emFields, layout);

            // get the particle velocity from t=n to t=n+1
            accelerate_(rangeOut, rangeOut, mass);

            // now advance the particles from t=n+1/2 to t=n+1 using v_{n+1} just calculated
            // and get a pointer to the first leaving particle
            firstLeaving = pushStep_(rangeOut, rangeOut, particleIsNotLeaving1, PushStep::PostPush);

            rangeOut = makeRange(rangeOut.begin(), std::move(firstLeaving));

            return rangeOut.end();
        }



        /** see Pusher::move() documentation*/
        virtual void setMeshAndTimeStep(std::array<double, dim> ms, double ts) override
        {
            std::transform(std::begin(ms), std::end(ms), std::begin(halfDtOverDl_),
                           [ts](double& x) { return 0.5 * ts / x; });
            dt_ = ts;
        }



    private:
        enum class PushStep { PrePush, PostPush };

        /** move the particle partIn of half a time step and store it in partOut
         */
        template<typename Particle, typename ParticleIter>
        void advancePosition_(Particle const& partIn, ParticleIter& partOut)
        {
            std::array<int, dim> outCell;
            // push the particle
            for (std::size_t iDim = 0; iDim < dim; ++iDim)
            {
                double delta = partIn.delta[iDim]
                               + static_cast<double>(halfDtOverDl_[iDim] * partIn.v[iDim]);

                double iCell = std::floor(delta);
                if (std::abs(delta) > 2)
                {
                    PHARE_LOG_ERROR("Error, particle moves more than 1 cell, delta >2");
                }
                partOut->delta[iDim] = delta - iCell;
                partOut->iCell[iDim] = static_cast<int>(iCell + partIn.iCell[iDim]);
            }
            // partOut.change_icell(outCell);
        }



        /** advance the particles in rangeIn of half a time step and store them
         * in rangeOut.
         * @return the function returns and iterator on the first leaving particle, as
         * detected by the ParticleSelector
         */
        template<typename ParticleRangeIn, typename ParticleRangeOut>
        void pushStep_(ParticleRangeIn const& rangeIn, ParticleRangeOut& rangeOut, PushStep step)
        {
            auto currentOut = rangeOut.begin();
            for (auto& currentIn : rangeIn)
            {
                // in the first push, this is the first time
                // we push to rangeOut, which contains crap
                // the push will only touch the particle position
                // but the next step being the acceleration of
                // rangeOut, we need to copy rangeIn weights, charge
                // and velocity. This is done here although
                // not strictly speaking this function's business
                // to take advantage that we're already looping
                // over rangeIn particles.
                if (step == PushStep::PrePush)
                {
                    currentOut->charge = currentIn.charge;
                    currentOut->weight = currentIn.weight;
                    currentOut->v      = currentIn.v;
                }
                advancePosition_(currentIn, currentOut);
                ++currentOut;
            }
        }

        template<typename ParticleRangeIn, typename ParticleRangeOut>
        auto pushStep_(ParticleRangeIn const& rangeIn, ParticleRangeOut& rangeOut,
                       ParticleSelector const& particleIsNotLeaving, PushStep step)
        {
            pushStep_(rangeIn, rangeOut, step);

            // now all particles have been pushed
            // those not satisfying the predicate after the push
            // are found in [newEnd:end[
            // those for which pred is true are in [firstOut,newEnd[
            return std::partition(std::begin(rangeOut), std::end(rangeOut), particleIsNotLeaving);
        }




        /** Accelerate the particles in rangeIn and put the new velocity in rangeOut
         */
        template<typename ParticleRangeIn, typename ParticleRangeOut>
        void accelerate_(ParticleRangeIn inputParticles, ParticleRangeOut outputParticles,
                         double mass)
        {
            double dto2m = 0.5 * dt_ / mass;

            auto currentOut = outputParticles.begin();

            for (auto const& currentIn : inputParticles)
            {
                double coef1 = currentIn.charge * dto2m;

                // We now apply the 3 steps of the BORIS PUSHER

                // 1st half push of the electric field
                double velx1 = currentIn.v[0] + coef1 * currentIn.Ex;
                double vely1 = currentIn.v[1] + coef1 * currentIn.Ey;
                double velz1 = currentIn.v[2] + coef1 * currentIn.Ez;


                // preparing variables for magnetic rotation
                double const rx = coef1 * currentIn.Bx;
                double const ry = coef1 * currentIn.By;
                double const rz = coef1 * currentIn.Bz;

                double const rx2  = rx * rx;
                double const ry2  = ry * ry;
                double const rz2  = rz * rz;
                double const rxry = rx * ry;
                double const rxrz = rx * rz;
                double const ryrz = ry * rz;

                double const invDet = 1. / (1. + rx2 + ry2 + rz2);

                // preparing rotation matrix due to the magnetic field
                // m = invDet*(I + r*r - r x I) - I where x denotes the cross product
                double const mxx = 1. + rx2 - ry2 - rz2;
                double const mxy = 2. * (rxry + rz);
                double const mxz = 2. * (rxrz - ry);

                double const myx = 2. * (rxry - rz);
                double const myy = 1. + ry2 - rx2 - rz2;
                double const myz = 2. * (ryrz + rx);

                double const mzx = 2. * (rxrz + ry);
                double const mzy = 2. * (ryrz - rx);
                double const mzz = 1. + rz2 - rx2 - ry2;

                // magnetic rotation
                double const velx2 = (mxx * velx1 + mxy * vely1 + mxz * velz1) * invDet;
                double const vely2 = (myx * velx1 + myy * vely1 + myz * velz1) * invDet;
                double const velz2 = (mzx * velx1 + mzy * vely1 + mzz * velz1) * invDet;


                // 2nd half push of the electric field
                velx1 = velx2 + coef1 * currentIn.Ex;
                vely1 = vely2 + coef1 * currentIn.Ey;
                velz1 = velz2 + coef1 * currentIn.Ez;

                // Update particle velocity
                currentOut->v[0] = velx1;
                currentOut->v[1] = vely1;
                currentOut->v[2] = velz1;

                ++currentOut;
            }
        }




        std::array<double, dim> halfDtOverDl_;
        double dt_;
    };

} // namespace core

} // namespace PHARE


#endif