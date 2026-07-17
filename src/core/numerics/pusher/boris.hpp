#ifndef PHARE_CORE_PUSHER_BORIS_HPP
#define PHARE_CORE_PUSHER_BORIS_HPP

#include "core/errors.hpp"
#include "core/logger.hpp"
#include "core/vector.hpp"
#include "core/utilities/range/range.hpp"
#include "core/data/particles/particle.hpp"
#include "core/data/particles/particle_array_def.hpp"


#include <array>
#include <cmath>
#include <cstddef>
#include <algorithm>
#include <iterator>
#include <stdexcept>

namespace PHARE::core
{


template<std::size_t dim, typename ParticleRange, typename Electromag, typename Interpolator,
         typename BoundaryCondition, typename GridLayout>
class BorisPusher
{
public:
    using ParticleSelector = std::function<ParticleRange(ParticleRange&)>;


    BorisPusher() {} // default for shared_ptr usage
    BorisPusher(std::array<double, dim> const& ms, double const ts) { setMeshAndTimeStep(ms, ts); }

    // This move function should be considered when being used so that all particles are pushed
    // twice - see: https://github.com/PHAREHUB/PHARE/issues/571
    /** see Pusher::move() documentation*/
#if 0
    ParticleRange move(ParticleRange const& rangeIn, ParticleRange& rangeOut,
                       Electromag const& emFields, double mass, Interpolator& interpolator,
                       ParticleSelector const& particleIsNotLeaving, BoundaryCondition& bc,
                       GridLayout const& layout)
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
#endif


    ParticleRange move(ParticleRange const& rangeIn, ParticleRange& rangeOut,
                       Electromag const& emFields, double mass, Interpolator& interpolator,
                       GridLayout const& layout, ParticleSelector firstSelector,
                       ParticleSelector secondSelector)
    {
        if (rangeIn.size() == 0)
            return rangeOut;

        PHARE_LOG_SCOPE(2, "Boris::move_no_bc");

        // push the particles of half a step
        // rangeIn : t=n, rangeOut : t=n+1/2
        // Do not partition on this step - this is to keep all domain and ghost
        //   particles consistent. see: https://github.com/PHAREHUB/PHARE/issues/571
        prePushStep_(rangeIn, rangeOut);

        rangeOut = firstSelector(rangeOut);

        double const dto2m = 0.5 * dt_ / mass;

        for (auto idx = rangeOut.ibegin(); idx < rangeOut.iend(); ++idx)
        {
            auto& particles = rangeOut.array();

            //  get electromagnetic fields interpolated on the particles of rangeOut stop at
            //  newEnd. get the particle velocity from t=n to t=n+1
            accelerate_(particles, interpolator(particles, emFields, layout, idx), dto2m, idx);

            // now advance the particles from t=n+1/2 to t=n+1 using v_{n+1} just calculated
            // and get a pointer to the first leaving particle
            postPushStep_(rangeOut.array(), idx, halfDtOverDl_);
        }

        return secondSelector(rangeOut);
    }



    /** see Pusher::move() documentation*/
    void setMeshAndTimeStep(std::array<double, dim> const& ms, double const ts) _PHARE_ALL_FN_
    {
        std::transform(std::begin(ms), std::end(ms), std::begin(halfDtOverDl_),
                       [ts](double const& x) { return 0.5 * ts / x; });
        dt_ = ts;
    }



private:
    /** move the particle partIn of half a time step and store it in partOut
     */
    template<typename ParticleIn, typename ParticleOut>
    auto static advancePosition_(ParticleIn const& partIn, ParticleOut& partOut,
                                 std::array<double, dim> halfDtOverDl) _PHARE_ALL_FN_
    {
        std::array<int, dim> newCell;
        for (std::size_t iDim = 0; iDim < dim; ++iDim)
        {
            double delta
                = partIn.delta()[iDim] + static_cast<double>(halfDtOverDl[iDim] * partIn.v()[iDim]);

            double iCell = std::floor(delta);
            if (std::abs(delta) > 2)
            {
                // PHARE_LOG_ERROR("Error, particle moves more than 1 cell, delta >2");
                throw_runtime_error("Error, particle moves more than 1 cell, delta >2");
            }
            partOut.delta()[iDim] = delta - iCell;
            newCell[iDim]         = static_cast<int>(iCell + partIn.iCell()[iDim]);
        }
        return newCell;
    }



    /** advance the particles in rangeIn of half a time step and store them
     * in rangeOut.
     * @return the function returns and iterator on the first leaving particle, as
     * detected by the ParticleSelector
     */
    void prePushStep_(ParticleRange const& rangeIn, ParticleRange& rangeOut)
    {
        using ParticleArray = std::decay_t<decltype(rangeIn.array())>;

        auto& inParticles  = rangeIn.array();
        auto& outParticles = rangeOut.array();
        for (auto inIdx = rangeIn.ibegin(), outIdx = rangeOut.ibegin(); inIdx < rangeIn.iend();
             ++inIdx, ++outIdx)
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

            outParticles.charge(outIdx) = inParticles.charge(inIdx);
            outParticles.weight(outIdx) = inParticles.weight(inIdx);
            outParticles.v(outIdx)      = inParticles.v(inIdx);

            auto in      = inParticles.begin() + inIdx;
            auto out     = outParticles.begin() + outIdx;
            auto newCell = advancePosition_(in, out, halfDtOverDl_);

            if constexpr (any_in(ParticleArray::layout_mode, LayoutMode::AoSMapped))
            {
                if (newCell != inParticles.iCell(inIdx))
                    outParticles.change_icell(newCell, outIdx);
            }
            else
                outParticles.iCell(outIdx) = newCell;
        }
    }

    template<typename Particles>
    void static postPushStep_(Particles& particles, std::size_t idx,
                              std::array<double, dim> halfDtOverDl) _PHARE_ALL_FN_
    {
        auto particle = particles.begin() + idx;
        auto newCell  = advancePosition_(deref(particle), deref(particle), halfDtOverDl);

        if constexpr (any_in(Particles::layout_mode, LayoutMode::AoSMapped))
        {
            if (newCell != particles.iCell(idx))
                particles.change_icell(newCell, idx);
        }
        else
            particles.iCell(idx) = newCell;
    }


    /** Accelerate the particles in rangeIn and put the new velocity in rangeOut
     */

    template<typename Particles, typename ParticleEB>
    void static accelerate_(Particles& particles, ParticleEB const& particleEB, double const& dto2m,
                            std::size_t const idx) _PHARE_ALL_FN_
    {
        auto& [pE, pB]        = particleEB;
        auto& [pEx, pEy, pEz] = pE;
        auto& [pBx, pBy, pBz] = pB;

        double coef1 = particles.charge(idx) * dto2m;

        // We now apply the 3 steps of the BORIS PUSHER

        // 1st half push of the electric field
        auto const& in_v                = particles.v(idx);
        std::array<double, 3> const vel = {in_v[0] + coef1 * pEx, //
                                           in_v[1] + coef1 * pEy, //
                                           in_v[2] + coef1 * pEz};

        // preparing variables for magnetic rotation
        double const rx = coef1 * pBx;
        double const ry = coef1 * pBy;
        double const rz = coef1 * pBz;

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
        double const velx2 = (mxx * vel[0] + mxy * vel[1] + mxz * vel[2]) * invDet;
        double const vely2 = (myx * vel[0] + myy * vel[1] + myz * vel[2]) * invDet;
        double const velz2 = (mzx * vel[0] + mzy * vel[1] + mzz * vel[2]) * invDet;

        // 2nd half push of the electric field
        // Update particle velocity
        particles.v(idx) = {velx2 + coef1 * pEx, //
                            vely2 + coef1 * pEy, //
                            velz2 + coef1 * pEz};
    }




    std::array<double, dim> halfDtOverDl_;
    double dt_;
};

} // namespace PHARE::core


#endif
