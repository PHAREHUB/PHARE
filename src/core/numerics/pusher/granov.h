#ifndef PHARE_CORE_PUSHER_GRANOV_H
#define PHARE_CORE_PUSHER_GRANOV_H


#include "core/ThreadPool.h"

#include <array>
#include <cmath>
#include <thread>
#include <cstddef>
#include <iostream>


#include "core/numerics/pusher/pusher.h"
#include "core/utilities/range/range.h"
#include "core/logger.h"
#include "core/def.h"

namespace PHARE::core
{
template<std::size_t dim, typename ParticleIterator, typename Electromag, typename Interpolator,
         typename BoundaryCondition, typename GridLayout>
class GranovPusher
    : public Pusher<dim, ParticleIterator, Electromag, Interpolator, BoundaryCondition, GridLayout>
{
public:
    using Super
        = Pusher<dim, ParticleIterator, Electromag, Interpolator, BoundaryCondition, GridLayout>;
    using ParticleSelector = typename Super::ParticleSelector;
    using ParticleRange    = Range<ParticleIterator>;

    GranovPusher() _PHARE_FN_SIG_ {}

    GranovPusher(GridLayout const& layout, double timestep, double mass) _PHARE_FN_SIG_
    {
        setMeshAndTimeStep(layout.meshSize(), timestep);
        accelerate_setup(mass);
    }
    ~GranovPusher() _PHARE_FN_SIG_ {}

    /** see Pusher::move() documentation*/
    ParticleIterator move(ParticleRange const& rangeIn, ParticleRange& rangeOut,
                          Electromag const& emFields, double mass, Interpolator& interpolator,
                          ParticleSelector const& particleIsNotLeaving, BoundaryCondition& bc,
                          GridLayout const& layout) override
    {
        throw_runtime_error("NOT TO BE CALLED");
        // return rangeOut.end();
    }


    /** see Pusher::move() documentation*/
    ParticleIterator move(ParticleRange const& rangeIn, ParticleRange& rangeOut,
                          Electromag const& emFields, double mass, Interpolator& interpolator,
                          ParticleSelector const& particleIsNotLeaving,
                          GridLayout const& layout) override
    {
        throw_runtime_error("NOT TO BE CALLED");
        // return rangeOut.end();
    }


    template<typename Particle_t, typename Selector>
    bool move_in_place(Particle_t& particle, Electromag const& emFields, Interpolator& interpolator,
                       Selector const& particleIsNotLeaving,
                       GridLayout const& layout) _PHARE_FN_SIG_
    {
        advancePosition_(particle);

        auto eb = interpolator.meshToParticle(particle, emFields, layout);

        accelerate_(particle, eb);

        advancePosition_(particle);

        return particleIsNotLeaving(particle);
    }


    /** see Pusher::move() documentation*/
    ParticleIterator move(ParticleRange& range, Electromag const& emFields, double mass,
                          Interpolator& interpolator, ParticleSelector const& particleIsNotLeaving,
                          GridLayout const& layout) override
    {
        return range.end();
    }




    /** see Pusher::move() documentation*/
    virtual void setMeshAndTimeStep(std::array<double, dim> const& ms,
                                    double ts) override _PHARE_FN_SIG_
    {
        // std::transform(std::begin(ms), std::end(ms), std::begin(halfDtOverDl_),
        //                [ts](double const& x) { return 0.5 * ts / x; });
        for (std::size_t i = 0; i < dim; ++i)
        {
            halfDtOverDl_[i] = 0.5 * ts / ms[i];
            assert(halfDtOverDl_[i] != 0);
        }
        // halfDtOverDl_ = ms;
        dt_ = ts;
        assert(dt_ > 0);
    }


    void accelerate_setup(double mass) _PHARE_FN_SIG_
    {
        dto2m_ = 0.5 * dt_ / mass;
        assert(dto2m_ > 0);
    }

private:
    enum class PushStep { PrePush, PostPush };

    /** move the particle partIn of half a time step and store it in partOut
     */
    template<typename ParticleIter>
    void advancePosition_(ParticleIter& particle) _PHARE_FN_SIG_
    {
        advancePosition_(particle, particle);
    }

    template<typename Particle>
    void advancePosition_(Particle const& partIn, Particle& partOut) _PHARE_FN_SIG_
    {
        // push the particle
        for (std::size_t iDim = 0; iDim < dim; ++iDim)
        {
            double delta = partIn.delta[iDim] + (halfDtOverDl_[iDim] * partIn.v[iDim]);

            double iCell = std::floor(delta);
            //             assert(!(std::abs(delta) > 2)); //
            assert(delta > -2);
            assert(delta < 2);

            partOut.delta[iDim] = delta - iCell;
            assert(partOut.iCell[iDim] > -5);
            partOut.iCell[iDim] = static_cast<int>(iCell + partIn.iCell[iDim]);
            assert(partOut.iCell[iDim] > -5);
        }
    }


    /** Accelerate the particles in rangeIn and put the new velocity in rangeOut
     */
    template<typename Particle_t, typename ParticleElectromag>
    void accelerate_(Particle_t& particle,
                     ParticleElectromag const& particleElectromag) _PHARE_FN_SIG_
    {
        accelerate_(particle, particle, particleElectromag);
    }

    template<typename Particle_t, typename ParticleElectromag>
    void accelerate_(Particle_t const& currentIn, Particle_t& currentOut,
                     ParticleElectromag const& particleElectromag) _PHARE_FN_SIG_
    {
        double coef1 = currentIn.charge * dto2m_;

        // We now apply the 3 steps of the BORIS PUSHER

        auto const& [Ex, Ey, Ez] = particleElectromag[0];
        auto const& [Bx, By, Bz] = particleElectromag[1];

        // 1st half push of the electric field
        double velx1 = currentIn.v[0] + coef1 * Ex;
        double vely1 = currentIn.v[1] + coef1 * Ey;
        double velz1 = currentIn.v[2] + coef1 * Ez;


        // preparing variables for magnetic rotation
        double const rx = coef1 * Bx;
        double const ry = coef1 * By;
        double const rz = coef1 * Bz;

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
        currentOut.v[0] = velx2 + coef1 * Ex;
        currentOut.v[1] = vely2 + coef1 * Ey;
        currentOut.v[2] = velz2 + coef1 * Ez;
    }




    std::array<double, dim> halfDtOverDl_{PHARE::core::ConstArray<double, dim>(1)};
    double dt_ = 1;

    double dto2m_ = 1; // not thread safe, cannot // multiple pops at the same time
};

} // namespace PHARE::core

#endif
