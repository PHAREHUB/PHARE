#ifndef PHARE_CORE_PUSHER_BORIS_HPP
#define PHARE_CORE_PUSHER_BORIS_HPP


#include "core/errors.hpp"
#include "core/logger.hpp"
#include "core/numerics/pusher/pusher.hpp"


#include <array>
#include <cmath>
#include <cstddef>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <string>
#include <type_traits>
#include <vector>

namespace PHARE::core
{


template<std::size_t dim, typename ParticleRange, typename Electromag, typename Interpolator,
         typename BoundaryCondition, typename GridLayout>
class BorisPusher
    : public Pusher<dim, ParticleRange, Electromag, Interpolator, BoundaryCondition, GridLayout>
{
public:
    using Super
        = Pusher<dim, ParticleRange, Electromag, Interpolator, BoundaryCondition, GridLayout>;

private:
    using ParticleSelector = Super::ParticleSelector;

public:
    // This move function should be considered when being used so that all particles are pushed
    // twice - see: https://github.com/PHAREHUB/PHARE/issues/571
    /** see Pusher::move() documentation*/
#if 0
    ParticleRange move(ParticleRange const& rangeIn, ParticleRange& rangeOut,
                       Electromag const& emFields, double charge, double mass,
                       Interpolator& interpolator,
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
#endif


    ParticleRange move(ParticleRange const& rangeIn, ParticleRange& rangeOut,
                       Electromag const& emFields, double charge, double mass,
                       Interpolator& interpolator, GridLayout const& layout,
                       ParticleSelector firstSelector, ParticleSelector secondSelector) override
    {
        PHARE_LOG_SCOPE(3, "Boris::move_no_bc");

        // push the particles of half a step
        // rangeIn : t=n, rangeOut : t=n+1/2
        // Do not partition on this step - this is to keep all domain and ghost
        //   particles consistent. see: https://github.com/PHAREHUB/PHARE/issues/571
        prePushStep_(rangeIn, rangeOut);

        rangeOut = firstSelector(rangeOut);

        double const coef1      = accelerationCoef(charge, mass);
        auto const halfDtOverDl = halfDtOverDl_;
        auto& particles         = rangeOut.array();
        auto const end          = rangeOut.iend();
        for (auto idx = rangeOut.ibegin(); idx < end; ++idx)
        {
            auto& currPart = particles[idx];

            //  get electromagnetic fields interpolated on the particles of rangeOut stop at newEnd.
            //  get the particle velocity from t=n to t=n+1
            auto const& local_em = interpolator(currPart, emFields, layout);
            accelerate_(currPart, local_em, coef1);

            // now advance the particles from t=n+1/2 to t=n+1 using v_{n+1} just calculated
            std::array<int, dim> newCell;
            if (!advancePosition_(halfDtOverDl, currPart, currPart, newCell)) [[unlikely]]
                throwMovedTwoCells_("PostPush", halfDtOverDl, currPart, &local_em, &layout);
            if (newCell != currPart.iCell)
                particles.change_icell(newCell, idx);
        }

        return secondSelector(rangeOut);
    }



    /** move() with a no-op first selector, pushing range in place in a single pass.
     * The cell map receives the same updates, in the same order, as move(): all half-step
     * crossings in index order, then all full-step crossings in index order.
     * Kept out of line: inlined into the caller, the compiler forms different FMAs than in
     * move(), and the result is no longer bitwise identical.
     */
    [[gnu::noinline]] ParticleRange moveInPlace(ParticleRange& range, Electromag const& emFields,
                                                double charge, double mass,
                                                Interpolator& interpolator,
                                                GridLayout const& layout,
                                                ParticleSelector const& selector)
    {
        PHARE_LOG_SCOPE(3, "Boris::moveInPlace");

        double const coef1      = accelerationCoef(charge, mass);
        auto const halfDtOverDl = halfDtOverDl_;
        auto& particles         = range.array();
        auto const end          = range.iend();

        preCrossings_.clear();
        postCrossings_.clear();

        for (auto idx = range.ibegin(); idx < end; ++idx)
        {
            auto& currPart = particles[idx];

            std::array<int, dim> midCell;
            if (!advancePosition_(halfDtOverDl, currPart, currPart, midCell)) [[unlikely]]
                throwMovedTwoCells_("PrePush", halfDtOverDl, currPart);
            if (midCell != currPart.iCell)
            {
                preCrossings_.push_back({idx, currPart.iCell, midCell});
                currPart.iCell = midCell;
            }

            auto const& local_em = interpolator(currPart, emFields, layout);
            accelerate_(currPart, local_em, coef1);

            std::array<int, dim> newCell;
            if (!advancePosition_(halfDtOverDl, currPart, currPart, newCell)) [[unlikely]]
                throwPostPushAfterPrePush_(range, idx, local_em, layout);
            if (newCell != midCell)
            {
                postCrossings_.push_back({idx, midCell, newCell});
                currPart.iCell = newCell;
            }
        }

        for (auto const& [idx, oldCell, newCell] : preCrossings_)
            particles.map_icell_change(oldCell, newCell, idx);
        for (auto const& [idx, oldCell, newCell] : postCrossings_)
            particles.map_icell_change(oldCell, newCell, idx);

        return selector(range);
    }



    /** charge * dt / (2 * mass), the coefficient of E and B in the Boris acceleration */
    double accelerationCoef(double const charge, double const mass) const
    {
        double const dto2m = 0.5 * dt_ / mass;
        return charge * dto2m;
    }


    /** push a single particle in place from t=n to t=n+1, with the same arithmetic as move()
     * but without selection or cell map update: the caller owns the particle (e.g. a local
     * copy). coef1 is accelerationCoef(charge, mass).
     */
    template<typename Particle_t>
    void push(Particle_t& part, Electromag const& emFields, Interpolator& interpolator,
              GridLayout const& layout, double const coef1) const
    {
        auto const halfDtOverDl = halfDtOverDl_;

        std::array<int, dim> newCell;
        if (!advancePosition_(halfDtOverDl, part, part, newCell)) [[unlikely]]
            throwMovedTwoCells_("PrePush", halfDtOverDl, part);
        part.iCell = newCell;

        auto const& local_em = interpolator(part, emFields, layout);
        accelerate_(part, local_em, coef1);

        if (!advancePosition_(halfDtOverDl, part, part, newCell)) [[unlikely]]
            throwMovedTwoCells_("PostPush", halfDtOverDl, part, &local_em, &layout);
        part.iCell = newCell;
    }


    /** see Pusher::move() documentation*/
    void setMeshAndTimeStep(std::array<double, dim> ms, double const ts) override
    {
        std::transform(std::begin(ms), std::end(ms), std::begin(halfDtOverDl_),
                       [ts](double& x) { return 0.5 * ts / x; });
        dt_ = ts;
    }



private:
    /** move the particle partIn of half a time step and store it in partOut.
     * Returns false, leaving partOut untouched, if it would cross more than 2 cells.
     */
    template<typename Particle>
    static bool advancePosition_(std::array<double, dim> const& halfDtOverDl,
                                 Particle const& partIn, Particle& partOut,
                                 std::array<int, dim>& newCell)
    {
        std::array<double, dim> delta;
        bool tooFar = false;
        for (std::size_t iDim = 0; iDim < dim; ++iDim)
        {
            delta[iDim]
                = partIn.delta[iDim] + static_cast<double>(halfDtOverDl[iDim] * partIn.v[iDim]);
            tooFar |= std::abs(delta[iDim]) > 2;
        }
        if (tooFar)
            return false;

        for (std::size_t iDim = 0; iDim < dim; ++iDim)
        {
            auto const iCell    = static_cast<int>(std::floor(delta[iDim]));
            partOut.delta[iDim] = delta[iDim] - iCell;
            newCell[iDim]       = iCell + partIn.iCell[iDim];
        }
        return true;
    }


    // error path of advancePosition_, kept out of the push loops
    template<typename Particle, typename ParticleEB = std::nullptr_t>
    [[noreturn, gnu::noinline, gnu::cold]] static void
    throwMovedTwoCells_(char const* step, std::array<double, dim> const& halfDtOverDl,
                        Particle const& part, ParticleEB const* particleEB = nullptr,
                        GridLayout const* layout = nullptr)
    {
        std::stringstream ss;
        for (std::size_t iDim = 0; iDim < dim; ++iDim)
        {
            double const delta
                = part.delta[iDim] + static_cast<double>(halfDtOverDl[iDim] * part.v[iDim]);
            if (std::abs(delta) > 2)
            {
                ss << step << " Particle moved 2 cells with delta/vel: ";
                ss << delta << "/" << part.v[iDim] << std::endl;
                break;
            }
        }

        DictionaryException ex{"cause", ss.str()};
        if constexpr (!std::is_same_v<ParticleEB, std::nullptr_t>)
        {
            auto const& [e, b] = *particleEB;
            for (std::uint16_t i = 0; i < 3; ++i)
                ex("E_" + std::to_string(i), std::to_string(e[i]));
            for (std::uint16_t i = 0; i < 3; ++i)
                ex("B_" + std::to_string(i), std::to_string(b[i]));
            ex("level", std::to_string(layout->levelNumber()));
        }
        throw ex;
    }


    /** moveInPlace error path: move() would report the first half-step failure of any particle
     * before a full-step failure, so check the remaining half steps first.
     */
    template<typename ParticleEB>
    [[noreturn, gnu::noinline, gnu::cold]] void
    throwPostPushAfterPrePush_(ParticleRange const& range, std::size_t const failedIdx,
                               ParticleEB const& particleEB, GridLayout const& layout) const
    {
        auto const em         = particleEB;
        auto const& particles = range.array();
        for (auto idx = failedIdx + 1; idx < range.iend(); ++idx)
        {
            auto part = particles[idx];
            std::array<int, dim> cell;
            if (!advancePosition_(halfDtOverDl_, part, part, cell))
                throwMovedTwoCells_("PrePush", halfDtOverDl_, part);
        }
        throwMovedTwoCells_("PostPush", halfDtOverDl_, particles[failedIdx], &em, &layout);
    }


    /** advance the particles in rangeIn of half a time step and store them
     * in rangeOut.
     */
    void prePushStep_(ParticleRange const& rangeIn, ParticleRange& rangeOut)
    {
        auto const halfDtOverDl = halfDtOverDl_;
        auto& inParticles       = rangeIn.array();
        auto& outParticles      = rangeOut.array();
        auto const inEnd        = rangeIn.iend();
        for (auto inIdx = rangeIn.ibegin(), outIdx = rangeOut.ibegin(); inIdx < inEnd;
             ++inIdx, ++outIdx)
        {
            // in the first push, this is the first time
            // we push to rangeOut, which contains crap
            // the push will only touch the particle position
            // but the next step being the acceleration of
            // rangeOut, we need to copy rangeIn weights
            // and velocity. This is done here although
            // not strictly speaking this function's business
            // to take advantage that we're already looping
            // over rangeIn particles.

            outParticles[outIdx].weight = inParticles[inIdx].weight;
            outParticles[outIdx].v      = inParticles[inIdx].v;

            std::array<int, dim> newCell;
            if (!advancePosition_(halfDtOverDl, inParticles[inIdx], outParticles[outIdx], newCell))
                [[unlikely]]
                throwMovedTwoCells_("PrePush", halfDtOverDl, inParticles[inIdx]);
            if (newCell != inParticles[inIdx].iCell)
                outParticles.change_icell(newCell, outIdx);
        }
    }


    /** Accelerate the particles in rangeIn and put the new velocity in rangeOut
     */
    template<typename Particle_t, typename ParticleEB>
    static void accelerate_(Particle_t& part, ParticleEB const& particleEB, double const coef1)
    {
        auto& [pE, pB]        = particleEB;
        auto& [pEx, pEy, pEz] = pE;
        auto& [pBx, pBy, pBz] = pB;

        // We now apply the 3 steps of the BORIS PUSHER

        // 1st half push of the electric field
        double velx1 = part.v[0] + coef1 * pEx;
        double vely1 = part.v[1] + coef1 * pEy;
        double velz1 = part.v[2] + coef1 * pEz;


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
        double const velx2 = (mxx * velx1 + mxy * vely1 + mxz * velz1) * invDet;
        double const vely2 = (myx * velx1 + myy * vely1 + myz * velz1) * invDet;
        double const velz2 = (mzx * velx1 + mzy * vely1 + mzz * velz1) * invDet;


        // 2nd half push of the electric field
        velx1 = velx2 + coef1 * pEx;
        vely1 = vely2 + coef1 * pEy;
        velz1 = velz2 + coef1 * pEz;

        // Update particle velocity
        part.v[0] = velx1;
        part.v[1] = vely1;
        part.v[2] = velz1;
    }




    std::array<double, dim> halfDtOverDl_;
    double dt_;

    // moveInPlace cell changes, replayed on the cell map after the push
    struct CellCrossing
    {
        std::size_t idx;
        std::array<int, dim> oldCell, newCell;
    };
    std::vector<CellCrossing> preCrossings_, postCrossings_;
};

} // namespace PHARE::core


#endif
