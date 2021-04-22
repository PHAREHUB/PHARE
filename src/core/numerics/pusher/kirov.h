#ifndef PHARE_CORE_PUSHER_KIROV_H
#define PHARE_CORE_PUSHER_KIROV_H


#include "core/ThreadPool.h"

#include <array>
#include <cmath>
#include <thread>
#include <cstddef>
#include <iostream>

#include "core/numerics/pusher/pusher.h"
#include "core/utilities/range/range.h"
#include "core/logger.h"

namespace PHARE::core
{
template<std::size_t dim, typename ParticleIterator, typename Electromag, typename Interpolator,
         typename BoundaryCondition, typename GridLayout, typename ThreadPool_t = EXT::ThreadPool>
class KirovPusher
    : public Pusher<dim, ParticleIterator, Electromag, Interpolator, BoundaryCondition, GridLayout>
{
public:
    using Super
        = Pusher<dim, ParticleIterator, Electromag, Interpolator, BoundaryCondition, GridLayout>;
    using ParticleSelector = typename Super::ParticleSelector;
    using ParticleRange    = Range<ParticleIterator>;

    KirovPusher(std::size_t threads = 1)
        : pool_{threads - 1}
    {
        assert(threads > 0);
    }

    /** see Pusher::move() documentation*/
    ParticleIterator move(ParticleRange const& rangeIn, ParticleRange& rangeOut,
                          Electromag const& emFields, double mass, Interpolator& interpolator,
                          ParticleSelector const& particleIsNotLeaving, BoundaryCondition& bc,
                          GridLayout const& layout) override
    {
        // push the particles of half a step
        // rangeIn : t=n, rangeOut : t=n+1/Z
        // get a pointer on the first particle of rangeOut that leaves the patch
        auto firstLeaving = pushStep_(rangeIn, rangeOut, particleIsNotLeaving, PushStep::PrePush);

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
        PHARE_LOG_SCOPE("Kirov::move_no_bc");

        // push the particles of half a step
        // rangeIn : t=n, rangeOut : t=n+1/2
        // get a pointer on the first particle of rangeOut that leaves the patch
        auto firstLeaving = pushStep_(rangeIn, rangeOut, particleIsNotLeaving, PushStep::PrePush);

        rangeOut = makeRange(rangeOut.begin(), std::move(firstLeaving));

        // get electromagnetic fields interpolated on the particles of rangeOut
        // stop at newEnd.
        interpolator(rangeOut.begin(), rangeOut.end(), emFields, layout);

        // get the particle velocity from t=n to t=n+1
        accelerate_(rangeOut, rangeOut, mass);

        // now advance the particles from t=n+1/2 to t=n+1 using v_{n+1} just calculated
        // and get a pointer to the first leaving particle
        firstLeaving = pushStep_(rangeOut, rangeOut, particleIsNotLeaving, PushStep::PostPush);

        rangeOut = makeRange(rangeOut.begin(), std::move(firstLeaving));

        return rangeOut.end();
    }


    template<typename Particle_t>
    bool move_in_place(Particle_t& particle, Electromag const& emFields, Interpolator& interpolator,
                       ParticleSelector const& particleIsNotLeaving, GridLayout const& layout)
    {
        advancePosition_(particle);

        if (particleIsNotLeaving(particle))
        {
            interpolator.meshToParticle(particle, emFields, layout);

            accelerate_(particle);

            advancePosition_(particle);

            return particleIsNotLeaving(particle);
        }
        return false;
    }


    /** see Pusher::move() documentation*/
    auto move(ParticleRange& range, Electromag const& emFields, double mass,
              Interpolator& interpolator, ParticleSelector const& particleIsNotLeaving,
              GridLayout const& layout) /*override*/
    {
        PHARE_LOG_SCOPE("Kirov::move_in_place_no_bc");

        accelerate_setup(mass);

        std::vector<bool> particlesAreNotLeaving(range.size(), false);

        auto _move = [&](std::size_t const from, std::size_t const to) {
            for (std::size_t i = from; i < to; ++i)
                particlesAreNotLeaving[i]
                    = move_in_place(range[i], emFields, interpolator, particleIsNotLeaving, layout);
        };

        std::size_t perThread = range.size() / (pool_.size() + 1);
        for (std::size_t i = 0; i < pool_.size(); ++i)
            pool_.enqueue([&, ix = i]() {
                auto from = ix * perThread;
                auto to   = from + perThread;
                _move(from, to);
            });

        _move(pool_.size() * perThread, range.size()); // should start at 0 for cache?

        pool_.sync();

        auto firstLeaving
            = std::partition(std::begin(range), std::end(range), [&](auto const& part) {
                  return particlesAreNotLeaving[&part - &(*range.begin())];
              });

        return makeRange(range.begin(), std::move(firstLeaving)).end();
    }




    /** see Pusher::move() documentation*/
    virtual void setMeshAndTimeStep(std::array<double, dim> ms, double ts) override
    {
        std::transform(std::begin(ms), std::end(ms), std::begin(halfDtOverDl_),
                       [ts](double const& x) { return 0.5 * ts / x; });
        dt_ = ts;
    }



private:
    enum class PushStep { PrePush, PostPush };

    /** move the particle partIn of half a time step and store it in partOut
     */
    template<typename ParticleIter>
    void advancePosition_(ParticleIter& particle)
    {
        advancePosition_(particle, particle);
    }

    template<typename ParticleIter>
    void advancePosition_(ParticleIter const& partIn, ParticleIter& partOut)
    {
        // push the particle
        for (std::size_t iDim = 0; iDim < dim; ++iDim)
        {
            float delta
                = partIn.delta[iDim] + static_cast<float>(halfDtOverDl_[iDim] * partIn.v[iDim]);

            float iCell = std::floor(delta);
            if (std::abs(delta) > 2)
            {
                throw std::runtime_error("Error, particle moves more than 1 cell, delta >2");
            }
            partOut.delta[iDim] = delta - iCell;
            partOut.iCell[iDim] = static_cast<int>(iCell + partIn.iCell[iDim]);
        }
    }



    /** advance the particles in rangeIn of half a time step and store them
     * in rangeOut.
     * @return the function returns and iterator on the first leaving particle, as
     * detected by the ParticleSelector
     */
    template<typename ParticleRangeIn, typename ParticleRangeOut>
    auto pushStep_(ParticleRangeIn const& rangeIn, ParticleRangeOut& rangeOut,
                   ParticleSelector const& particleIsNotLeaving, PushStep step)
    {
        auto currentOut = rangeOut.begin();

        for (auto& currentIn : rangeIn)
        {
            if (step == PushStep::PrePush)
            {
                currentOut->charge = currentIn.charge;
                currentOut->weight = currentIn.weight;
                currentOut->v      = currentIn.v;
            }
            // push the particle
            advancePosition_(currentIn, *currentOut);
            currentOut++;
        }

        // now all particles have been pushed
        // those not satisfying the predicate after the push
        // are found in [newEnd:end[
        // those for which pred is true are in [firstOut,newEnd[
        return std::partition(std::begin(rangeOut), std::end(rangeOut), particleIsNotLeaving);
    }



    void accelerate_setup(double mass) { dto2m_ = 0.5 * dt_ / mass; }

    template<typename ParticleRangeIn, typename ParticleRangeOut>
    void accelerate_(ParticleRangeIn inputParticles, ParticleRangeOut outputParticles, double mass)
    {
        accelerate_setup(mass);

        auto currentOut = outputParticles.begin();

        for (auto const& currentIn : inputParticles)
        {
            accelerate_(currentIn, *currentOut);
            ++currentOut;
        }
    }


    /** Accelerate the particles in rangeIn and put the new velocity in rangeOut
     */
    template<typename Particle_t>
    void accelerate_(Particle_t& particle)
    {
        accelerate_(particle, particle);
    }

    template<typename Particle_t>
    void accelerate_(Particle_t const& currentIn, Particle_t& currentOut)
    {
        double coef1 = currentIn.charge * dto2m_;

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

        // 2nd half push of the electric field / Update particle velocity
        currentOut.v[0] = velx2 + coef1 * currentIn.Ex;
        currentOut.v[1] = vely2 + coef1 * currentIn.Ey;
        currentOut.v[2] = velz2 + coef1 * currentIn.Ez;
    }




    std::array<double, dim> halfDtOverDl_;
    double dt_;

    double dto2m_; // not thread safe, cannot // multiple pops at the same time
    ThreadPool_t pool_;
};

} // namespace PHARE::core

#endif
