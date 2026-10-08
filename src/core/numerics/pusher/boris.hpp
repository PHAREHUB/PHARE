#ifndef PHARE_CORE_PUSHER_BORIS_HPP
#define PHARE_CORE_PUSHER_BORIS_HPP

#include "core/logger.hpp"
#include "core/errors.hpp"
#include "core/numerics/pusher/boris/basics.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include <array>
#include <cstddef>
#include <sstream>
#include <iterator>
#include <algorithm>

namespace PHARE::core
{

template<std::size_t dim>
class BorisPusher
{
public:
    BorisPusher() {} // default for shared_ptr usage
    BorisPusher(std::array<double, dim> const& ms, double const ts) { setMeshAndTimeStep(ms, ts); }

    // This move function should be considered when being used so that all particles are pushed
    // twice - see: https://github.com/PHAREHUB/PHARE/issues/571
    /** see Pusher::move() documentation*/
#if 0
    auto move(auto const& rangeIn, auto& rangeOut,
                       auto const& emFields, double mass, auto& interpolator,
                       auto const& particleIsNotLeaving, BoundaryCondition& bc,
                       auto const& layout)
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


    auto move(auto const& rangeIn, auto& rangeOut, auto const& emFields, double mass,
              auto& interpolator, auto const& layout, auto firstSelector, auto secondSelector)
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

            auto const local_em = interpolator(particles, emFields, layout, idx);
            accelerate_(particles, local_em, dto2m, idx);

            try
            {
                postPushStep_(rangeOut.array(), idx, halfDtOverDl_);
            }
            catch (DictionaryException const& bex)
            {
                auto ex             = bex;
                auto const& [e, b] = local_em;
                for (std::uint16_t i = 0; i < 3; ++i)
                    ex("E_" + std::to_string(i), std::to_string(e[i]));
                for (std::uint16_t i = 0; i < 3; ++i)
                    ex("B_" + std::to_string(i), std::to_string(b[i]));
                ex("level", std::to_string(layout.levelNumber()));
                throw ex;
            }
        }

        return secondSelector(rangeOut);
    }


    /** see Pusher::move() documentation*/
    void setMeshAndTimeStep(std::array<double, dim> const& ms, double const ts)
    {
        std::transform(std::begin(ms), std::end(ms), std::begin(halfDtOverDl_),
                       [ts](double const& x) { return 0.5 * ts / x; });
        dt_ = ts;
    }


private:
    /** advance the particles in rangeIn of half a time step and store them
     * in rangeOut.
     * @return the function returns and iterator on the first leaving particle, as
     * detected by the auto
     */
    void prePushStep_(auto const& rangeIn, auto& rangeOut)
    {
        using ParticleArray              = std::decay_t<decltype(rangeIn.array())>;
        static constexpr auto alloc_mode = ParticleArray::alloc_mode;

        auto& inParticles  = rangeIn.array();
        auto& outParticles = rangeOut.array();
        for (auto inIdx = rangeIn.ibegin(), outIdx = rangeOut.ibegin(); inIdx < rangeIn.iend();
             ++inIdx, ++outIdx)
        {
            outParticles.charge(outIdx) = inParticles.charge(inIdx);
            outParticles.weight(outIdx) = inParticles.weight(inIdx);
            outParticles.v(outIdx)      = inParticles.v(inIdx);
            outParticles.delta(outIdx)  = inParticles.delta(inIdx);
            outParticles.iCell(outIdx)  = inParticles.iCell(inIdx);

            auto out = outParticles.begin() + outIdx;

            std::array<int, dim> newCell;
            try
            {
                newCell = boris::advance<alloc_mode>(deref(out), halfDtOverDl_);
            }
            catch (boris::MoveTwoCellException const& e)
            {
                std::stringstream ss;
                ss << "PrePush Particle moved 2 cells with delta/vel: ";
                ss << e.delta << "/" << e.vel;
                throw DictionaryException{}("cause", ss.str());
            }

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
                              std::array<double, dim> halfDtOverDl)
    {
        static constexpr auto alloc_mode = Particles::alloc_mode;
        auto particle                    = particles.begin() + idx;

        std::array<int, dim> newCell;
        try
        {
            newCell = boris::advance<alloc_mode>(deref(particle), halfDtOverDl);
        }
        catch (boris::MoveTwoCellException const& e)
        {
            std::stringstream ss;
            ss << "PostPush Particle moved 2 cells with delta/vel: ";
            ss << e.delta << "/" << e.vel;
            throw DictionaryException{}("cause", ss.str());
        }

        if constexpr (any_in(Particles::layout_mode, LayoutMode::AoSMapped))
        {
            if (newCell != particles.iCell(idx))
                particles.change_icell(newCell, idx);
        }
        else
            particles.iCell(idx) = newCell;
    }

    template<typename Particles, typename ParticleEB>
    void static accelerate_(Particles& particles, ParticleEB const& particleEB, double const& dto2m,
                            std::size_t const idx)
    {
        auto particle = particles.begin() + idx;
        boris::accelerate(deref(particle), particleEB, dto2m);
    }


    std::array<double, dim> halfDtOverDl_;
    double dt_;
};

} // namespace PHARE::core


#endif
