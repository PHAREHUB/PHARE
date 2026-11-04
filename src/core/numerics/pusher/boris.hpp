#ifndef PHARE_CORE_PUSHER_BORIS_HPP
#define PHARE_CORE_PUSHER_BORIS_HPP


#include "core/errors.hpp"
#include "core/logger.hpp"
#include "core/numerics/pusher/pusher.hpp"
#include "core/numerics/pusher/particle_lanes.hpp"
#include "core/numerics/pusher/hyb_xsimd.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/models/quantities/hybrid_quantities.hpp"


#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <sstream>
#include <stdexcept>
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
    ParticleRange moveInPlace(ParticleRange& range, Electromag const& emFields, double charge,
                              double mass, Interpolator& interpolator, GridLayout const& layout,
                              ParticleSelector const& selector)
    {
        PHARE_LOG_SCOPE(3, "Boris::moveInPlace");
        return moveInPlace_<hyb::Backend>(range, emFields, charge, mass, interpolator, layout,
                                          selector);
    }

    template<typename Backend>
    [[gnu::noinline]] ParticleRange
    moveInPlace_(ParticleRange& range, Electromag const& emFields, double charge, double mass,
                 Interpolator& interpolator, GridLayout const& layout,
                 ParticleSelector const& selector)
    {
        double const coef1      = accelerationCoef(charge, mass);
        auto const halfDtOverDl = halfDtOverDl_;
        auto& particles         = range.array();
        auto const end          = range.iend();

        preCrossings_.clear();
        postCrossings_.clear();

        auto moveOne = [&](std::size_t const idx) {
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
        };

        auto const fields = blockFields_<Interpolator>(emFields, layout);
        PushBlock block;
        // xsimd: crossings and lanes written back by writeBackLanesAoS_ (register transpose),
        // where it is faster than the per-lane write-back: GCC in 2D and 3D, clang in 3D
#if PHARE_HYB_XSIMD and defined(__clang__)
        constexpr bool vecWB = dim == 3;
#elif PHARE_HYB_XSIMD
        constexpr bool vecWB = dim > 1;
#else
        constexpr bool vecWB = false;
#endif
        pushBlocks_<Interpolator, Backend>(
            particles, range.ibegin(), end, fields, coef1, block,
            [&](std::size_t const idx, std::size_t const n, bool const ok) {
                if (!ok) [[unlikely]]
                {
                    // a particle of the block moves more than one cell: the scalar path throws
                    // the same exception as move()
                    for (auto i = idx; i < idx + n; ++i)
                        moveOne(i);
                    return;
                }
#if PHARE_HYB_XSIMD
                if constexpr (vecWB)
                {
                    writeBackLanesAoS_(particles, idx, n, block);
                    return;
                }
#endif
                for (std::size_t l = 0; l < n; ++l)
                {
                    auto& part = particles[idx + l];
                    std::array<int, dim> oldCell, midCell, newCell;
                    for (std::size_t iDim = 0; iDim < dim; ++iDim)
                    {
                        oldCell[iDim]    = part.iCell[iDim];
                        midCell[iDim]    = block.midCell[iDim][l];
                        newCell[iDim]    = block.cell[iDim][l];
                        part.delta[iDim] = block.delta[iDim][l];
                    }
                    for (std::size_t c = 0; c < 3; ++c)
                        part.v[c] = block.v[c][l];
                    if (midCell != oldCell)
                        preCrossings_.push_back({idx + l, oldCell, midCell});
                    if (newCell != midCell)
                        postCrossings_.push_back({idx + l, midCell, newCell});
                    part.iCell = newCell;
                }
            });

        if (particles.mapValid()) // a sorted array does not keep its cell map
        {
            for (auto const& [idx, oldCell, newCell] : preCrossings_)
                particles.map_icell_change(oldCell, newCell, idx);
            for (auto const& [idx, oldCell, newCell] : postCrossings_)
                particles.map_icell_change(oldCell, newCell, idx);
        }

        return selector(range);
    }


    /** fn(idx, cellBefore, cellAfter) for every particle of the last moveInPlace whose iCell
     * differs between before and after the push, in ascending idx (the input of the bag sort,
     * ParticleArray::sortBags)
     */
    template<typename Fn>
    void forEachCellChange(Fn&& fn) const
    {
        std::size_t ipre = 0, ipost = 0;
        while (ipre < preCrossings_.size() or ipost < postCrossings_.size())
        {
            bool const havePre  = ipre < preCrossings_.size();
            bool const havePost = ipost < postCrossings_.size();
            if (havePre and havePost and preCrossings_[ipre].idx == postCrossings_[ipost].idx)
            {
                if (preCrossings_[ipre].oldCell != postCrossings_[ipost].newCell)
                    fn(preCrossings_[ipre].idx, preCrossings_[ipre].oldCell,
                       postCrossings_[ipost].newCell);
                ++ipre, ++ipost;
            }
            else if (havePre
                     and (!havePost or preCrossings_[ipre].idx < postCrossings_[ipost].idx))
            {
                fn(preCrossings_[ipre].idx, preCrossings_[ipre].oldCell,
                   preCrossings_[ipre].newCell);
                ++ipre;
            }
            else
            {
                fn(postCrossings_[ipost].idx, postCrossings_[ipost].oldCell,
                   postCrossings_[ipost].newCell);
                ++ipost;
            }
        }
    }



    /** use the cell kernels (same-cell gather, pushBlockCell_) on the sorted arrays pushed
     * next; set by the ion updater per domain array
     */
    void cellKernels(bool const on) { cellKernels_ = on; }

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


    /** push() on a copy of every particle of particles, in order, calling onPushed(pushed)
     * after each. Particles are pushed by blocks of blockSize with the same per-particle
     * arithmetic as push(); a block where a particle moves more than one cell is redone with
     * push(), which throws.
     */
    template<typename Particles, typename OnPushed>
    void pushEach(Particles const& particles, Electromag const& emFields,
                  Interpolator& interpolator, GridLayout const& layout, double const coef1,
                  OnPushed&& onPushed)
    {
        pushEach_<hyb::Backend>(particles, emFields, interpolator, layout, coef1, onPushed);
    }

    template<typename Backend, typename Particles, typename OnPushed>
    [[gnu::noinline]] void pushEach_(Particles const& particles, Electromag const& emFields,
                                     Interpolator& interpolator, GridLayout const& layout,
                                     double const coef1, OnPushed&& onPushed)
    {
        using Particle_t = std::decay_t<decltype(particles[0])>;

        auto pushOne = [&](std::size_t const idx) {
            Particle_t pushed = particles[idx];
            push(pushed, emFields, interpolator, layout, coef1);
            onPushed(pushed);
        };

        auto const fields = blockFields_<Interpolator>(emFields, layout);
        auto const end    = particles.size();
        PushBlock block;
        pushBlocks_<Interpolator, Backend>(
            particles, 0, end, fields, coef1, block,
            [&](std::size_t const idx, std::size_t const n, bool const ok) {
                if (!ok) [[unlikely]]
                {
                    for (auto i = idx; i < idx + n; ++i)
                        pushOne(i);
                    return;
                }
                for (std::size_t l = 0; l < n; ++l)
                {
                    Particle_t pushed = particles[idx + l];
                    for (std::size_t iDim = 0; iDim < dim; ++iDim)
                    {
                        pushed.iCell[iDim] = block.cell[iDim][l];
                        pushed.delta[iDim] = block.delta[iDim][l];
                    }
                    for (std::size_t c = 0; c < 3; ++c)
                        pushed.v[c] = block.v[c][l];
                    onPushed(pushed);
                }
            });
    }


    /** pushEach with the deposit of the pushed copies done by blocks (interp order >= 2):
     * a pushed block is deposited straight from its lanes, lane l only if keep(new AMR cell),
     * by Interpolator::depositLanesOOL, in particle order. targets come from
     * Interpolator::depositTargets.
     */
    template<typename Particles, typename Targets, typename Keep>
    void pushEachDeposit(Particles const& particles, Electromag const& emFields,
                         Interpolator& interpolator, GridLayout const& layout, double const coef1,
                         Targets const& targets, Keep const& keep)
    {
        pushEachDeposit_<hyb::Backend>(particles, emFields, interpolator, layout, coef1, targets,
                                       keep);
    }

#if PHARE_HYB_XSIMD
    /** pushEach with the deposit of the pushed copies through a CellDeposit window (sorted
     * arrays, cell kernels): blocks per cell row, lanes masked by keep(new AMR cell), the
     * window keyed on the row's cell (pre-push)
     */
    template<typename Particles, typename Targets, typename CellDep, typename Keep>
    void pushEachDepositCell(Particles const& particles, Electromag const& emFields,
                             Interpolator& interpolator, GridLayout const& layout,
                             double const coef1, Targets const& targets, CellDep& cellDep,
                             Keep const& keep)
    {
        pushEachDepositCell_<hyb::Backend>(particles, emFields, interpolator, layout, coef1,
                                           targets, cellDep, keep);
    }

    template<typename Backend, typename Particles, typename Targets, typename CellDep,
             typename Keep>
    [[gnu::noinline]] void
    pushEachDepositCell_(Particles const& particles, Electromag const& emFields,
                         Interpolator& interpolator, GridLayout const& layout, double const coef1,
                         Targets const& targets, CellDep& cellDep, Keep const& keep)
    {
        using Particle_t = std::decay_t<decltype(particles[0])>;

        auto const fields = blockFields_<Interpolator>(emFields, layout);
        auto const end    = particles.size();
        PushBlock block;
        double weight[blockSize];
        int mask[blockSize];
        hyb::LanePtrs<dim> lanes;
        for (std::size_t iDim = 0; iDim < dim; ++iDim)
        {
            lanes.delta[iDim] = block.delta[iDim];
            lanes.cell[iDim]  = block.cell[iDim];
        }
        for (std::size_t c = 0; c < 3; ++c)
            lanes.v[c] = block.v[c];
        lanes.weight = weight;

        cellDep.begin(targets, fields.amrToLocal, /*withMargin=*/true);
        // blocks per cell row (runs = true): the window is flushed once per cell, and the
        // block's first particle (pre-push) gives the row's cell
        pushBlocks_<Interpolator, Backend>(
            particles, 0, end, fields, coef1, block,
            [&](std::size_t const idx, std::size_t const n, bool const ok) {
                if (!ok) [[unlikely]]
                {
                    for (auto i = idx; i < idx + n; ++i)
                    {
                        Particle_t pushed = particles[i];
                        push(pushed, emFields, interpolator, layout, coef1);
                        if (keep(pushed.iCell))
                            interpolator.depositOn(pushed, targets, layout);
                    }
                    return;
                }
                std::array<int, dim> key;
                for (std::size_t iDim = 0; iDim < dim; ++iDim)
                    key[iDim] = particles[idx].iCell[iDim];
                for (std::size_t l = 0; l < blockSize; ++l)
                {
                    weight[l] = particles[idx + std::min(l, n - 1)].weight;
                    std::array<int, dim> cell;
                    for (std::size_t iDim = 0; iDim < dim; ++iDim)
                        cell[iDim] = block.cell[iDim][l];
                    mask[l] = l < n and keep(cell);
                }
                cellDep.block(lanes, mask, key);
            },
            /*runs=*/true);
        cellDep.end();
    }
#endif // PHARE_HYB_XSIMD

    template<typename Backend, typename Particles, typename Targets, typename Keep>
    [[gnu::noinline]] void pushEachDeposit_(Particles const& particles, Electromag const& emFields,
                                            Interpolator& interpolator, GridLayout const& layout,
                                            double const coef1, Targets const& targets,
                                            Keep const& keep)
    {
        using Particle_t = std::decay_t<decltype(particles[0])>;

        auto const fields = blockFields_<Interpolator>(emFields, layout);
        auto const end    = particles.size();
        PushBlock block;
        double weight[blockSize];
        hyb::LanePtrs<dim> lanes;
        for (std::size_t iDim = 0; iDim < dim; ++iDim)
        {
            lanes.delta[iDim] = block.delta[iDim];
            lanes.cell[iDim]  = block.cell[iDim];
        }
        for (std::size_t c = 0; c < 3; ++c)
            lanes.v[c] = block.v[c];
        lanes.weight = weight;

        pushBlocks_<Interpolator, Backend>(
            particles, 0, end, fields, coef1, block,
            [&](std::size_t const idx, std::size_t const n, bool const ok) {
                if (!ok) [[unlikely]]
                {
                    for (auto i = idx; i < idx + n; ++i)
                    {
                        Particle_t pushed = particles[i];
                        push(pushed, emFields, interpolator, layout, coef1);
                        if (keep(pushed.iCell))
                            interpolator.depositOn(pushed, targets, layout);
                    }
                    return;
                }
                for (std::size_t l = 0; l < blockSize; ++l)
                    weight[l] = particles[idx + std::min(l, n - 1)].weight;
                Interpolator::template depositLanesOOL<Backend>(lanes, n, fields.amrToLocal,
                                                                targets, 1., keep);
            });
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
        accelerate_(part.v[0], part.v[1], part.v[2], pEx, pEy, pEz, pBx, pBy, pBz, coef1);
    }

    [[gnu::always_inline]] static void accelerate_(double& vx, double& vy, double& vz,
                                                   double const pEx, double const pEy,
                                                   double const pEz, double const pBx,
                                                   double const pBy, double const pBz,
                                                   double const coef1)
    {
        // We now apply the 3 steps of the BORIS PUSHER

        // the coef1 * E products feed both half pushes; with implicit contraction, the compiler
        // fuses them or the rotation products depending on the inlining context (vectorised or
        // not), which changes the result by round-off. Explicit fma pins one choice.

        // 1st half push of the electric field
        double velx1 = std::fma(coef1, pEx, vx);
        double vely1 = std::fma(coef1, pEy, vy);
        double velz1 = std::fma(coef1, pEz, vz);


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
        velx1 = std::fma(coef1, pEx, velx2);
        vely1 = std::fma(coef1, pEy, vely2);
        velz1 = std::fma(coef1, pEz, velz2);

        // Update particle velocity
        vx = velx1;
        vy = vely1;
        vz = velz1;
    }



    // ---- block push: blockSize particles at a time, lane loops the compiler vectorises ----

    static constexpr std::size_t blockSize = 8;

    // block state after the push, one lane per particle
    struct PushBlock
    {
        double delta[dim][blockSize];
        int midCell[dim][blockSize];
        int cell[dim][blockSize];
        double v[3][blockSize];
    };

    // raw storage of an E or B component: data and row strides of its local index
    struct BlockField
    {
        double const* data;
        std::array<std::uint32_t, dim> shape;
    };

    struct BlockFields
    {
        std::array<BlockField, 6> EB; // Ex, Ey, Ez, Bx, By, Bz
        std::array<int, dim> amrToLocal;
    };

    template<typename Interpolator_t>
    static BlockFields blockFields_(Electromag const& emFields, GridLayout const& layout)
    {
        BlockFields f;
        auto const& [Ex, Ey, Ez] = emFields.E();
        auto const& [Bx, By, Bz] = emFields.B();
        auto set                 = [](auto const& field, BlockField& bf) {
            bf.data = field.data();
            for (std::size_t iDim = 0; iDim < dim; ++iDim)
                bf.shape[iDim] = field.shape()[iDim];
        };
        set(Ex, f.EB[0]);
        set(Ey, f.EB[1]);
        set(Ez, f.EB[2]);
        set(Bx, f.EB[3]);
        set(By, f.EB[4]);
        set(Bz, f.EB[5]);
        // GridLayout::AMRToLocal for cells: local = AMR - (AMRBox.lower - localStart)
        int const localStart = layout.physicalStartIndex(QtyCentering::dual, Direction::X);
        for (std::size_t iDim = 0; iDim < dim; ++iDim)
            f.amrToLocal[iDim] = localStart - layout.AMRBox().lower[iDim];
        return f;
    }


    // half-step position advance of the block, lane by lane as advancePosition_. floor(delta)
    // is converted to int by comparisons (|delta| <= 2 on valid lanes), which vectorises without
    // -fno-trapping-math. Returns a nonzero value if a lane moves more than two cells.
    // Delta/Cell/V: lane arrays [dim][blockSize] or arrays of lane pointers
    template<typename Delta, typename Cell, typename V>
    [[gnu::always_inline]] static bool advanceBlock_(std::array<double, dim> const& halfDtOverDl,
                                                     Delta& delta, Cell& cell, V const& v)
    {
        // the lane loops stay loops (unroll 1) for the loop vectoriser: once fully unrolled,
        // the comparisons are not SLP vectorised
        int far[blockSize] = {};
        for (std::size_t iDim = 0; iDim < dim; ++iDim)
        {
            double const hdl = halfDtOverDl[iDim];
            double d[blockSize];
            int shift[blockSize];
#pragma GCC unroll 1
            for (std::size_t l = 0; l < blockSize; ++l)
                d[l] = delta[iDim][l] + static_cast<double>(hdl * v[iDim][l]);
#pragma GCC unroll 1
            for (std::size_t l = 0; l < blockSize; ++l)
            {
                far[l] |= std::abs(d[l]) > 2.;
                // floor(d) for |d| <= 2
                shift[l] = (d[l] >= 1.) + (d[l] >= 2.) - (d[l] < 0.) - (d[l] < -1.);
            }
            // ivdep: lane pointers are not known to be disjoint
#if !defined(__clang__) // clang has no GCC ivdep (-Wunknown-pragmas)
#pragma GCC ivdep
#endif
#pragma GCC unroll 1
            for (std::size_t l = 0; l < blockSize; ++l)
            {
                delta[iDim][l] = d[l] - static_cast<double>(shift[l]);
                cell[iDim][l] += shift[l];
            }
        }
        int anyFar = 0;
        for (std::size_t l = 0; l < blockSize; ++l)
            anyFar |= far[l];
        return anyFar != 0;
    }


    // interpolation of one component on the block lanes, same sum order as MeshToParticle
    template<std::size_t nbrPts, typename Starts, typename Weights>
    [[gnu::always_inline]] static void gatherBlock_(BlockField const& f, Starts const& starts,
                                                    Weights const& weights,
                                                    double (&out)[blockSize])
    {
        double const* __restrict__ data = f.data;
        if constexpr (dim == 1)
        {
            auto const& [xs] = starts;
            auto const& [xw] = weights;
            for (std::size_t l = 0; l < blockSize; ++l)
            {
                double acc = 0.;
                for (std::size_t ix = 0; ix < nbrPts; ++ix)
                    acc += data[(*xs)[l] + static_cast<int>(ix)] * (*xw)[ix][l];
                out[l] = acc;
            }
        }
        else if constexpr (dim == 2)
        {
            auto const& [xs, ys] = starts;
            auto const& [xw, yw] = weights;
            int const ny         = f.shape[1];
            for (std::size_t l = 0; l < blockSize; ++l)
            {
                double acc = 0.;
                for (std::size_t ix = 0; ix < nbrPts; ++ix)
                {
                    int const row = ((*xs)[l] + static_cast<int>(ix)) * ny + (*ys)[l];
                    double Y      = 0.;
                    for (std::size_t iy = 0; iy < nbrPts; ++iy)
                        Y += data[row + static_cast<int>(iy)] * (*yw)[iy][l];
                    acc += Y * (*xw)[ix][l];
                }
                out[l] = acc;
            }
        }
        else
        {
            auto const& [xs, ys, zs] = starts;
            auto const& [xw, yw, zw] = weights;
            int const ny             = f.shape[1];
            int const nz             = f.shape[2];
            for (std::size_t l = 0; l < blockSize; ++l)
            {
                double acc = 0.;
                for (std::size_t ix = 0; ix < nbrPts; ++ix)
                {
                    double Y = 0.;
                    for (std::size_t iy = 0; iy < nbrPts; ++iy)
                    {
                        int const row = (((*xs)[l] + static_cast<int>(ix)) * ny + (*ys)[l]
                                         + static_cast<int>(iy))
                                            * nz
                                        + (*zs)[l];
                        double Z = 0.;
                        for (std::size_t iz = 0; iz < nbrPts; ++iz)
                            Z += data[row + static_cast<int>(iz)] * (*zw)[iz][l];
                        Y += Z * (*yw)[iy][l];
                    }
                    acc += Y * (*xw)[ix][l];
                }
                out[l] = acc;
            }
        }
    }


    // start indices and weights of the block lanes, as Interpolator::localIndexAndWeights_
    template<typename Interpolator_t, QtyCentering centering, std::size_t nbrPts, typename Delta>
    [[gnu::always_inline]] static void
    indexAndWeightsBlock_(int const (&localCell)[dim][blockSize], Delta const& delta,
                          int (&starts)[dim][blockSize], double (&weights)[dim][nbrPts][blockSize])
    {
        constexpr auto order = Interpolator_t::interp_order;
        for (std::size_t iDim = 0; iDim < dim; ++iDim)
            for (std::size_t l = 0; l < blockSize; ++l)
            {
                double const d = delta[iDim][l];
                int const lc   = localCell[iDim][l];
                int const s
                    = lc
                      - Interpolator_t::template computeStartLeftShift<QtyCentering, centering>(d);
                double normalizedPos = lc + d;
                if constexpr (centering == QtyCentering::dual)
                    normalizedPos -= .5;
                std::array<double, nbrPts> w;
                Weighter<order>{}.computeWeight(normalizedPos, s, w);
                starts[iDim][l] = s;
                for (std::size_t k = 0; k < nbrPts; ++k)
                    weights[iDim][k][l] = w[k];
            }
    }


    /** push particles [first, first + n), n <= blockSize, from t=n to t=n+1 into the first n
     * lanes of block, with the arithmetic of push(); particles are not modified. Returns false if
     * a particle moves more than one cell in a half step.
     */
    // out of line: one copy of the block kernel, shared by the push loops
    template<typename Interpolator_t, typename Backend, typename Particles>
    [[gnu::noinline]] bool pushBlock_(Particles const& particles, std::size_t const first,
                                      std::size_t const n, BlockFields const& f, double const coef1,
                                      PushBlock& b) const
    {
        // a partial block (n < blockSize) repeats its last particle in the unused lanes, so that
        // every particle goes through the same, vectorised, arithmetic
        for (std::size_t l = 0; l < blockSize; ++l)
        {
            auto const& part = particles[first + std::min(l, n - 1)];
            for (std::size_t iDim = 0; iDim < dim; ++iDim)
            {
                b.delta[iDim][l] = part.delta[iDim];
                b.cell[iDim][l]  = part.iCell[iDim];
            }
            for (std::size_t c = 0; c < 3; ++c)
                b.v[c][l] = part.v[c];
        }

        return pushLanes_<Interpolator_t, Backend>(b, b.midCell, f, coef1);
    }


    /** push lanes [0, blockSize) of b (members delta, cell, v: lane arrays or lane pointers)
     * IN PLACE from t=n to t=n+1, with the arithmetic of push(); the cells after the first half
     * step go to midCell. Returns false if a lane moves more than one cell in a half step (b is
     * then partially updated).
     * Backend: vectorisation backend (hyb::AutoVec: lane loops left to the compiler;
     * hyb::Xsimd: pushLanesXs_).
     */
    template<typename Interpolator_t, typename Backend, typename Lanes, typename MidCell>
    [[gnu::always_inline]] inline bool pushLanes_(Lanes& b, MidCell& midCell, BlockFields const& f,
                                                  double const coef1) const
    {
#if PHARE_HYB_XSIMD
        if constexpr (std::is_same_v<Backend, hyb::Xsimd>)
            return pushLanesXs_<Interpolator_t>(b, midCell, f, coef1);
#endif
        static_assert(std::is_same_v<Backend, hyb::AutoVec> or std::is_same_v<Backend, hyb::Xsimd>);
        constexpr std::size_t nbrPts = nbrPointsSupport(Interpolator_t::interp_order);
        using Scalar                 = HybridQuantity::Scalar;

        auto const halfDtOverDl = halfDtOverDl_;

        bool tooFar = advanceBlock_(halfDtOverDl, b.delta, b.cell, b.v);

        int localCell[dim][blockSize];
        for (std::size_t iDim = 0; iDim < dim; ++iDim)
        {
            int const amrToLocal = f.amrToLocal[iDim]; // midCell may alias f for the compiler
#if !defined(__clang__)
#pragma GCC ivdep
#endif
            for (std::size_t l = 0; l < blockSize; ++l)
            {
                midCell[iDim][l]   = b.cell[iDim][l];
                localCell[iDim][l] = b.cell[iDim][l] + amrToLocal;
            }
        }

        int dStarts[dim][blockSize], pStarts[dim][blockSize];
        double dWeights[dim][nbrPts][blockSize], pWeights[dim][nbrPts][blockSize];
        indexAndWeightsBlock_<Interpolator_t, QtyCentering::dual>(localCell, b.delta, dStarts,
                                                                  dWeights);
        indexAndWeightsBlock_<Interpolator_t, QtyCentering::primal>(localCell, b.delta, pStarts,
                                                                    pWeights);

        double EB[6][blockSize];
        auto gather = [&]<auto quantity>(std::size_t const c) {
            constexpr auto centerings = GridLayout::centering(quantity);
            auto starts               = [&]<std::size_t... i>(std::index_sequence<i...>) {
                return std::make_tuple(
                    (centerings[i] == QtyCentering::primal ? &pStarts[i] : &dStarts[i])...);
            }(std::make_index_sequence<dim>{});
            auto weights = [&]<std::size_t... i>(std::index_sequence<i...>) {
                return std::make_tuple(
                    (centerings[i] == QtyCentering::primal ? &pWeights[i] : &dWeights[i])...);
            }(std::make_index_sequence<dim>{});
            gatherBlock_<nbrPts>(f.EB[c], starts, weights, EB[c]);
        };
        gather.template operator()<Scalar::Ex>(0);
        gather.template operator()<Scalar::Ey>(1);
        gather.template operator()<Scalar::Ez>(2);
        gather.template operator()<Scalar::Bx>(3);
        gather.template operator()<Scalar::By>(4);
        gather.template operator()<Scalar::Bz>(5);

#if !defined(__clang__)
#pragma GCC ivdep
#endif
        for (std::size_t l = 0; l < blockSize; ++l)
            accelerate_(b.v[0][l], b.v[1][l], b.v[2][l], EB[0][l], EB[1][l], EB[2][l], EB[3][l],
                        EB[4][l], EB[5][l], coef1);

        tooFar |= advanceBlock_(halfDtOverDl, b.delta, b.cell, b.v);
        return !tooFar;
    }



#if PHARE_HYB_XSIMD
    // ---- hyb::Xsimd backend: the lane core of pushLanes_ on explicit xsimd batches ----
    // same per-lane arithmetic and expression order as the AutoVec lane loops

    // half-step advance of the lanes (advanceBlock_), registers in and out
    [[gnu::always_inline]] static bool advanceXs_(std::array<double, dim> const& halfDtOverDl,
                                                  hyb::xs::bd (&delta)[dim],
                                                  hyb::xs::bi (&cell)[dim],
                                                  hyb::xs::bd const (&v)[3])
    {
        using namespace hyb::xs;
        bm far{false};
        for (std::size_t iDim = 0; iDim < dim; ++iDim)
        {
            bd const d = xsimd::fma(bd(halfDtOverDl[iDim]), v[iDim], delta[iDim]);
            far        = far | (xsimd::abs(d) > bd(2.));
            // floor(d) for |d| <= 2, by comparisons (exact in double)
            bd const shift
                = ones(d >= bd(1.)) + ones(d >= bd(2.)) - ones(d < bd(0.)) - ones(d < bd(-1.));
            delta[iDim] = d - shift;
            cell[iDim]  = cell[iDim] + toI(shift);
        }
        return xsimd::any(far);
    }

    /** one component interpolated on the lanes (gatherBlock_): the stencil rows along the last
     * (contiguous) index loaded as one row of 4 doubles per lane and transposed
     * (hyb::xs::loadRows), the fma chain as GCC contracts gatherBlock_ (bitwise with the AutoVec
     * backend). Row loads are faster than indexed loads (vgatherdpd) at all orders.
     */
    template<std::size_t nbrPts>
    [[gnu::always_inline]] static hyb::xs::bd
    gatherXs_(BlockField const& f, hyb::xs::bi const* (&starts)[dim],
              hyb::xs::bd const (*(&weights)[dim])[nbrPts])
    {
        constexpr int W = 4;
        using namespace hyb::xs;
        double const* data = f.data;
        int size           = 1;
        for (std::size_t i = 0; i < dim; ++i)
            size *= static_cast<int>(f.shape[i]);
        bi const last(static_cast<int>(nbrPts) - 1);
        bd acc(0.);
        bd cols[nbrPts];
        if constexpr (dim == 1)
        {
            bi const row = *starts[0];
            loadRows(data, row, size, xsimd::any(row > bi(size - W)), cols);
            for (std::size_t ix = 0; ix < nbrPts; ++ix)
                acc = xsimd::fma(cols[ix], (*weights[0])[ix], acc);
        }
        else if constexpr (dim == 2)
        {
            bi const ny(static_cast<int>(f.shape[1]));
            // rows grow with ix: the last one bounds them all
            bool const edge = xsimd::any((*starts[0] + last) * ny + *starts[1] > bi(size - W));
            for (std::size_t ix = 0; ix < nbrPts; ++ix)
            {
                bi const row = (*starts[0] + bi(static_cast<int>(ix))) * ny + *starts[1];
                loadRows(data, row, size, edge, cols);
                bd Y(0.);
                for (std::size_t iy = 0; iy < nbrPts; ++iy)
                    Y = xsimd::fma(cols[iy], (*weights[1])[iy], Y);
                acc = xsimd::fma(Y, (*weights[0])[ix], acc);
            }
        }
        else
        {
            bi const ny(static_cast<int>(f.shape[1]));
            bi const nz(static_cast<int>(f.shape[2]));
            bool const edge = xsimd::any(
                ((*starts[0] + last) * ny + *starts[1] + last) * nz + *starts[2] > bi(size - W));
            for (std::size_t ix = 0; ix < nbrPts; ++ix)
            {
                bd Y(0.);
                for (std::size_t iy = 0; iy < nbrPts; ++iy)
                {
                    bi const row = ((*starts[0] + bi(static_cast<int>(ix))) * ny + *starts[1]
                                    + bi(static_cast<int>(iy)))
                                       * nz
                                   + *starts[2];
                    loadRows(data, row, size, edge, cols);
                    bd Z(0.);
                    for (std::size_t iz = 0; iz < nbrPts; ++iz)
                        Z = xsimd::fma(cols[iz], (*weights[2])[iz], Z);
                    Y = xsimd::fma(Z, (*weights[1])[iy], Y);
                }
                acc = xsimd::fma(Y, (*weights[0])[ix], acc);
            }
        }
        return acc;
    }

    // Boris rotation of the lanes (accelerate_), explicit fma where accelerate_ has std::fma
    [[gnu::always_inline]] static void accelerateXs_(hyb::xs::bd (&v)[3],
                                                     hyb::xs::bd const (&EB)[6], double const coef)
    {
        using hyb::xs::bd;
        bd const coef1(coef);
        bd const one(1.), two(2.);

        bd velx1 = xsimd::fma(coef1, EB[0], v[0]);
        bd vely1 = xsimd::fma(coef1, EB[1], v[1]);
        bd velz1 = xsimd::fma(coef1, EB[2], v[2]);

        bd const rx = coef1 * EB[3];
        bd const ry = coef1 * EB[4];
        bd const rz = coef1 * EB[5];

        // the contractions GCC forms in the AutoVec lane loop (read off its code), made
        // explicit: every r_i r_j product fused into its sum, rotation rows as
        // fma(m_z, v_z, fma(m_x, v_x, m_y * v_y))
        using xsimd::fma;
        using xsimd::fms;
        using xsimd::fnma;

        bd const invDet = one / fma(rz, rz, fma(ry, ry, fma(rx, rx, one)));

        bd const mxx = fnma(rz, rz, fnma(ry, ry, fma(rx, rx, one)));
        bd const mxy = two * fma(rx, ry, rz);
        bd const mxz = two * fms(rx, rz, ry);

        bd const myx = two * fms(rx, ry, rz);
        bd const myy = fnma(rz, rz, fnma(rx, rx, fma(ry, ry, one)));
        bd const myz = two * fma(ry, rz, rx);

        bd const mzx = two * fma(rx, rz, ry);
        bd const mzy = two * fms(ry, rz, rx);
        bd const mzz = fnma(ry, ry, fnma(rx, rx, fma(rz, rz, one)));

        bd const velx2 = fma(mxz, velz1, fma(mxx, velx1, mxy * vely1)) * invDet;
        bd const vely2 = fma(myz, velz1, fma(myx, velx1, myy * vely1)) * invDet;
        bd const velz2 = fma(mzz, velz1, fma(mzx, velx1, mzy * vely1)) * invDet;

        v[0] = xsimd::fma(coef1, EB[0], velx2);
        v[1] = xsimd::fma(coef1, EB[1], vely2);
        v[2] = xsimd::fma(coef1, EB[2], velz2);
    }

    /** the xsimd lane core on registers: delta, cell, v pushed in place from t=n to t=n+1, cells
     * after the first half step to midCell. Returns false if a lane moves too far.
     * unionGather (cell kernels, pushBlockCell_): called at mid step with the lanes' cells,
     * deltas and weights; returns true if it interpolated EB itself, else the row gather runs.
     * Without it (nullptr) the code is that of the other kernels.
     */
    template<typename Interpolator_t, typename MidCell, typename UnionGather = std::nullptr_t>
    [[gnu::always_inline]] inline bool coreXs_(hyb::xs::bd (&delta)[dim], hyb::xs::bi (&cell)[dim],
                                               hyb::xs::bd (&v)[3], MidCell& midCell,
                                               BlockFields const& f, double const coef1,
                                               UnionGather&& unionGather = nullptr) const
    {
        using namespace hyb::xs;
        constexpr std::size_t nbrPts = nbrPointsSupport(Interpolator_t::interp_order);
        constexpr int order          = Interpolator_t::interp_order;
        using Scalar                 = HybridQuantity::Scalar;

        auto const halfDtOverDl = halfDtOverDl_;

        bool tooFar = advanceXs_(halfDtOverDl, delta, cell, v);

        bi dS[dim], pS[dim];
        bd dW[dim][nbrPts], pW[dim][nbrPts];
        for (std::size_t iDim = 0; iDim < dim; ++iDim)
        {
            sti(&midCell[iDim][0], cell[iDim]);
            bi const localCell = cell[iDim] + bi(f.amrToLocal[iDim]);
            startAndWeights<order, true>(localCell, delta[iDim], dS[iDim], dW[iDim]);
            startAndWeights<order, false>(localCell, delta[iDim], pS[iDim], pW[iDim]);
        }

        bd EB[6];
        // always_inline: GCC otherwise outlines the lambda
        auto gather = [&]<auto quantity>(std::size_t const c) __attribute__((always_inline)) {
            constexpr auto centerings = GridLayout::centering(quantity);
            bi const* starts[dim];
            bd const(*weights[dim])[nbrPts];
            for (std::size_t i = 0; i < dim; ++i)
            {
                bool const primal = centerings[i] == QtyCentering::primal;
                starts[i]         = primal ? &pS[i] : &dS[i];
                weights[i]        = primal ? &pW[i] : &dW[i];
            }
            EB[c] = gatherXs_<nbrPts>(f.EB[c], starts, weights);
        };
        bool unionDone = false;
        if constexpr (!std::is_same_v<std::decay_t<UnionGather>, std::nullptr_t>)
            unionDone = unionGather(cell, delta, dW, pW, EB);
        if (!unionDone)
        {
            gather.template operator()<Scalar::Ex>(0);
            gather.template operator()<Scalar::Ey>(1);
            gather.template operator()<Scalar::Ez>(2);
            gather.template operator()<Scalar::Bx>(3);
            gather.template operator()<Scalar::By>(4);
            gather.template operator()<Scalar::Bz>(5);
        }

        accelerateXs_(v, EB, coef1);

        tooFar |= advanceXs_(halfDtOverDl, delta, cell, v);
        return !tooFar;
    }

    /** pushLanes_ on xsimd batches: lanes loaded once from b, kept in registers through both
     * half steps, stored back once (b untouched if a lane moves too far)
     */
    template<typename Interpolator_t, typename Lanes, typename MidCell>
    [[gnu::always_inline]] inline bool pushLanesXs_(Lanes& b, MidCell& midCell,
                                                    BlockFields const& f, double const coef1) const
    {
        using namespace hyb::xs;
        bd delta[dim], v[3];
        bi cell[dim];
        for (std::size_t iDim = 0; iDim < dim; ++iDim)
        {
            delta[iDim] = ld(&b.delta[iDim][0]);
            cell[iDim]  = ldi(&b.cell[iDim][0]);
        }
        for (std::size_t c = 0; c < 3; ++c)
            v[c] = ld(&b.v[c][0]);

        if (!coreXs_<Interpolator_t>(delta, cell, v, midCell, f, coef1)) [[unlikely]]
            return false;

        for (std::size_t iDim = 0; iDim < dim; ++iDim)
        {
            st(&b.delta[iDim][0], delta[iDim]);
            sti(&b.cell[iDim][0], cell[iDim]);
        }
        for (std::size_t c = 0; c < 3; ++c)
            st(&b.v[c][0], v[c]);
        return true;
    }

    // ---- cell kernels: same-cell gather on sorted domain arrays ----
    //
    // Blocks are formed per cell row (pushBlocks_). At mid step, a block whose lanes all sit in
    // the row's cell is "uniform": per direction, the stencil start of a component is either
    // fixed or one of two values (lc - o - 1 + (delta >= .5)):
    //   order 1: primal fixed, dual shifts; order 2: primal shifts, dual fixed;
    //   order 3: primal fixed, dual shifts.
    // The per-direction union of the two stencils (nbrPts + 1 nodes where it shifts) is loaded
    // from the field once per cell and component (CellCache), and every lane picks its stencil
    // out of it by a select on delta < .5 along the shifting directions. The weights, the
    // products and the order of the fma chain are those of gatherXs_, only where the values come
    // from changes, so the result is bitwise that of gatherXs_.
    // The sort keys stay cells (not half cells): the rows, the box queries of the messengers and
    // the deposit windows are per cell, and half-cell keys would cut the lane fill of a row by
    // 2^dim. The select costs one blend per shifting direction, applied to the partial sums of
    // the outer directions (one chain per union row) and to the values of the last direction.

    // lc - (stencil start when delta < .5): the union starts there
    static constexpr int unionOffset_(int const order, bool const dual)
    {
        if (order == 1)
            return dual ? 1 : 0;
        if (order == 2)
            return 1;
        return dual ? 2 : 1;
    }

    template<int order, typename Centerings>
    static constexpr std::array<bool, dim> unionShift_(Centerings const& centerings)
    {
        std::array<bool, dim> shift{};
        for (std::size_t i = 0; i < dim; ++i)
            shift[i] = (centerings[i] == QtyCentering::dual) != (order == 2);
        return shift;
    }

    // the six field unions of one cell (union of the stencils of the particles of the cell)
    template<std::size_t nbrPts>
    struct CellCache
    {
        static constexpr std::size_t maxU = [] {
            std::size_t n = 1;
            for (std::size_t i = 0; i < dim; ++i)
                n *= nbrPts + 1;
            return n;
        }();
        alignas(64) double U[6][maxU];
        bool loaded = false;
    };

    // the union of the component on cell (local index lc): U[i0][i1][i2] row-major, sizes
    // nbrPts + Shift[d]
    template<std::size_t nbrPts, std::array<bool, dim> Shift, int order, typename Centerings>
    [[gnu::always_inline]] static void loadUnion_(BlockField const& f, Centerings const& centerings,
                                                  std::array<int, dim> const& lc, double* U)
    {
        std::array<int, dim> b;
        for (std::size_t i = 0; i < dim; ++i)
            b[i] = lc[i] - unionOffset_(order, centerings[i] == QtyCentering::dual);
        constexpr int U0   = nbrPts + Shift[0];
        double const* data = f.data;
        if constexpr (dim == 1)
        {
            for (int i = 0; i < U0; ++i)
                U[i] = data[b[0] + i];
        }
        else if constexpr (dim == 2)
        {
            constexpr int U1 = nbrPts + Shift[1];
            int const ny     = f.shape[1];
            for (int i = 0; i < U0; ++i)
                for (int j = 0; j < U1; ++j)
                    U[i * U1 + j] = data[(b[0] + i) * ny + b[1] + j];
        }
        else
        {
            constexpr int U1 = nbrPts + Shift[1];
            constexpr int U2 = nbrPts + Shift[2];
            int const ny     = f.shape[1];
            int const nz     = f.shape[2];
            for (int i = 0; i < U0; ++i)
                for (int j = 0; j < U1; ++j)
                    for (int k = 0; k < U2; ++k)
                        U[(i * U1 + j) * U2 + k]
                            = data[((b[0] + i) * ny + b[1] + j) * nz + b[2] + k];
        }
    }

    /** gatherXs_ from a union: chain of direction D over the union block U (directions >= D).
     * Same fma chain as gatherXs_: innermost acc = fma(value, w, acc) from 0, outer
     * acc = fma(partial sum, w, acc) from 0. On a shifting direction lane l uses union row
     * i + 1 - (delta_l < .5): for the last direction a select of the two values, for the outer
     * ones a select of the two partial sums (each computed with the same chain as gatherXs_).
     */
    template<std::size_t nbrPts, std::array<bool, dim> Shift, std::size_t D = 0>
    [[gnu::always_inline]] static hyb::xs::bd
    unionChainXs_(double const* U, hyb::xs::bm const (&m)[dim],
                  hyb::xs::bd const (*(&weights)[dim])[nbrPts])
    {
        using namespace hyb::xs;
        constexpr std::size_t stride = [] {
            std::size_t s = 1;
            for (std::size_t e = D + 1; e < dim; ++e)
                s *= nbrPts + Shift[e];
            return s;
        }();
        bd acc(0.);
        if constexpr (D == dim - 1)
        {
            for (std::size_t i = 0; i < nbrPts; ++i)
            {
                bd const val = Shift[D] ? xsimd::select(m[D], bd(U[i]), bd(U[i + 1])) : bd(U[i]);
                acc          = xsimd::fma(val, (*weights[D])[i], acc);
            }
        }
        else
        {
            constexpr std::size_t R = nbrPts + Shift[D];
            bd S[R];
            for (std::size_t r = 0; r < R; ++r)
                S[r] = unionChainXs_<nbrPts, Shift, D + 1>(U + r * stride, m, weights);
            for (std::size_t i = 0; i < nbrPts; ++i)
            {
                bd const val = Shift[D] ? xsimd::select(m[D], S[i], S[i + 1]) : S[i];
                acc          = xsimd::fma(val, (*weights[D])[i], acc);
            }
        }
        return acc;
    }

    /** the union gather of pushBlockCell_: if every lane's mid-step cell is rowCell, EB
     * interpolated from the cell's stencil unions (loaded into cache at the first uniform block
     * of the row) and true; else false (the caller's row gather runs)
     */
    template<typename Interpolator_t, std::size_t nbrPts, typename Cache>
    [[gnu::always_inline]] static bool
    unionGatherXs_(std::array<int, dim> const& rowCell, Cache& cache, BlockFields const& f,
                   hyb::xs::bi const (&midCell)[dim], hyb::xs::bd const (&midDelta)[dim],
                   hyb::xs::bd const (&dW)[dim][nbrPts], hyb::xs::bd const (&pW)[dim][nbrPts],
                   hyb::xs::bd (&EB)[6])
    {
        using namespace hyb::xs;
        constexpr int order = Interpolator_t::interp_order;
        using Scalar        = HybridQuantity::Scalar;
        auto same           = midCell[0] == bi(rowCell[0]);
        for (std::size_t i = 1; i < dim; ++i)
            same = same & (midCell[i] == bi(rowCell[i]));
        if (!xsimd::all(same))
            return false;

        auto forEachComponent = [&](auto&& fn) __attribute__((always_inline)) {
            fn.template operator()<Scalar::Ex>(0);
            fn.template operator()<Scalar::Ey>(1);
            fn.template operator()<Scalar::Ez>(2);
            fn.template operator()<Scalar::Bx>(3);
            fn.template operator()<Scalar::By>(4);
            fn.template operator()<Scalar::Bz>(5);
        };
        if (!cache.loaded) [[unlikely]]
        {
            std::array<int, dim> lc;
            for (std::size_t i = 0; i < dim; ++i)
                lc[i] = rowCell[i] + f.amrToLocal[i];
            forEachComponent([&]<auto quantity>(std::size_t const c) {
                constexpr auto centerings = GridLayout::centering(quantity);
                constexpr auto shift      = unionShift_<order>(centerings);
                loadUnion_<nbrPts, shift, order>(f.EB[c], centerings, lc, cache.U[c]);
            });
            cache.loaded = true;
        }

        bm m[dim];
        for (std::size_t i = 0; i < dim; ++i)
            m[i] = midDelta[i] < bd(.5);
        forEachComponent([&]<auto quantity>(std::size_t const c) __attribute__((always_inline)) {
            constexpr auto centerings = GridLayout::centering(quantity);
            constexpr auto shift      = unionShift_<order>(centerings);
            bd const(*weights[dim])[nbrPts];
            for (std::size_t i = 0; i < dim; ++i)
                weights[i] = centerings[i] == QtyCentering::primal ? &pW[i] : &dW[i];
            EB[c] = unionChainXs_<nbrPts, shift>(cache.U[c], m, weights);
        });
        return true;
    }

    /** the lanes of particles [first, first + n) (lane l >= n: particle n - 1) straight from the
     * AoS into registers: per lane one zmm load of the particle past its weight (iCell, delta,
     * v), an 8 x 8 transpose in registers, cells split out of their int pairs. Same values as
     * staging through a PushBlock, without the narrow stores and wide reloads (which cannot be
     * store-forwarded).
     */
    template<typename Particles>
    [[gnu::always_inline]] static void loadLanesAoS_(Particles const& particles,
                                                     std::size_t const first, std::size_t const n,
                                                     hyb::xs::bd (&delta)[dim],
                                                     hyb::xs::bi (&cell)[dim], hyb::xs::bd (&v)[3])
    {
        using namespace hyb::xs;
        using P                  = std::decay_t<decltype(particles[0])>;
        constexpr std::size_t qC = (dim + 1) / 2; // qwords of iCell (+ padding)
        constexpr std::size_t nq = qC + dim + 3;
        static_assert(offsetof(P, weight) == 0 and offsetof(P, iCell) == 8
                      and offsetof(P, delta) == 8 * (1 + qC)
                      and offsetof(P, v) == 8 * (1 + qC + dim) and sizeof(P) == 8 * (1 + nq));
        constexpr __mmask8 mask = static_cast<__mmask8>((1u << nq) - 1u);
        __m512d z[blockSize];
        // a full zmm read of a particle ends inside the next one: masked only for the last
        // particle of the array
        if (mask == 0xFF or first + n < particles.size()) [[likely]]
            for (std::size_t l = 0; l < blockSize; ++l)
                z[l] = _mm512_loadu_pd(
                    reinterpret_cast<double const*>(&particles[first + std::min(l, n - 1)]) + 1);
        else
            for (std::size_t l = 0; l < blockSize; ++l)
                z[l] = _mm512_maskz_loadu_pd(
                    mask,
                    reinterpret_cast<double const*>(&particles[first + std::min(l, n - 1)]) + 1);
        __m512d t[blockSize], u[blockSize];
        for (std::size_t p = 0; p < blockSize / 2; ++p)
        {
            t[2 * p]     = _mm512_unpacklo_pd(z[2 * p], z[2 * p + 1]);
            t[2 * p + 1] = _mm512_unpackhi_pd(z[2 * p], z[2 * p + 1]);
        }
        // u[4h + 0..3]: rows 4h..4h+3, columns (0|4), (2|6), (1|5), (3|7)
        for (std::size_t h = 0; h < 2; ++h)
        {
            u[4 * h + 0] = _mm512_shuffle_f64x2(t[4 * h], t[4 * h + 2], 0x88);
            u[4 * h + 1] = _mm512_shuffle_f64x2(t[4 * h], t[4 * h + 2], 0xDD);
            u[4 * h + 2] = _mm512_shuffle_f64x2(t[4 * h + 1], t[4 * h + 3], 0x88);
            u[4 * h + 3] = _mm512_shuffle_f64x2(t[4 * h + 1], t[4 * h + 3], 0xDD);
        }
        auto col = [&](std::size_t const q) __attribute__((always_inline)) {
            constexpr int src[8] = {0, 2, 1, 3, 0, 2, 1, 3}; // u index of column q (mod 4)
            return q < 4 ? _mm512_shuffle_f64x2(u[src[q]], u[4 + src[q]], 0x88)
                         : _mm512_shuffle_f64x2(u[src[q]], u[4 + src[q]], 0xDD);
        };
        for (std::size_t k = 0; k < qC; ++k)
        {
            __m512i const pair = _mm512_castpd_si512(col(k));
            cell[2 * k]        = bi(_mm512_cvtepi64_epi32(pair));
            if (2 * k + 1 < dim)
                cell[2 * k + 1] = bi(_mm512_cvtepi64_epi32(_mm512_srli_epi64(pair, 32)));
        }
        for (std::size_t iDim = 0; iDim < dim; ++iDim)
            delta[iDim] = bd(col(qC + iDim));
        for (std::size_t c = 0; c < 3; ++c)
            v[c] = bd(col(qC + dim + c));
    }

    /** pushBlock_ (xsimd) for a block of the row of rowCell: lanes loaded from the particles
     * (loadLanesAoS_); at mid step a uniform block interpolates from the cell's unions (loaded at
     * the first uniform block of the row), any other block takes the row gather of gatherXs_.
     */
    template<typename Interpolator_t, typename Particles, typename Cache>
    [[gnu::noinline]] bool pushBlockCell_(Particles const& particles, std::size_t const first,
                                          std::size_t const n,
                                          std::array<int, dim> const& rowCell, Cache& cache,
                                          BlockFields const& f, double const coef1,
                                          PushBlock& b) const
    {
        using namespace hyb::xs;
        constexpr std::size_t nbrPts = nbrPointsSupport(Interpolator_t::interp_order);

        bd delta[dim], v[3];
        bi cell[dim];
        loadLanesAoS_(particles, first, n, delta, cell, v);

        auto unionGather = [&](bi const (&midCell)[dim], bd const (&midDelta)[dim],
                               bd const (&dW)[dim][nbrPts], bd const (&pW)[dim][nbrPts],
                               bd (&EB)[6]) __attribute__((always_inline)) -> bool {
            return unionGatherXs_<Interpolator_t, nbrPts>(rowCell, cache, f, midCell, midDelta,
                                                          dW, pW, EB);
        };

        if (!coreXs_<Interpolator_t>(delta, cell, v, b.midCell, f, coef1, unionGather))
            [[unlikely]]
            return false;

        for (std::size_t iDim = 0; iDim < dim; ++iDim)
        {
            st(&b.delta[iDim][0], delta[iDim]);
            sti(&b.cell[iDim][0], cell[iDim]);
        }
        for (std::size_t c = 0; c < 3; ++c)
            st(&b.v[c][0], v[c]);
        return true;
    }

    /** corrector write-back of a pushed block (xsimd): crossings of lanes l < n recorded as
     * moveInPlace_ does (cells compared as int pairs: mid vs new in registers, old read from
     * the particle), then iCell, delta, v of particles [first, first + n) from an 8 x 8
     * register transpose of the lanes into particle rows, per particle one zmm store
     * (3D: iCell + its padding, delta, v = the 64 B past the weight) or a ymm + xmm store (2D:
     * the 48 B past the weight). The weight and the next particle are not touched.
     */
    template<typename Particles>
    [[gnu::always_inline]] void writeBackLanesAoS_(Particles& particles, std::size_t const first,
                                                   std::size_t const n, PushBlock const& b)
    {
        using P                  = std::decay_t<decltype(particles[0])>;
        constexpr std::size_t qC = (dim + 1) / 2;
        constexpr std::size_t nq = qC + dim + 3;
        static_assert(offsetof(P, weight) == 0 and offsetof(P, iCell) == 8
                      and offsetof(P, delta) == 8 * (1 + qC)
                      and offsetof(P, v) == 8 * (1 + qC + dim) and sizeof(P) == 8 * (1 + nq));
        static_assert(nq == 8 or nq == 6);
        __m512d x[blockSize];
        for (std::size_t k = 0; k < qC; ++k)
        {
            __m512i const lo = _mm512_cvtepu32_epi64(
                _mm256_loadu_si256(reinterpret_cast<__m256i const*>(b.cell[2 * k])));
            __m512i hi = _mm512_setzero_si512();
            if (2 * k + 1 < dim)
                hi = _mm512_slli_epi64(_mm512_cvtepu32_epi64(_mm256_loadu_si256(
                                           reinterpret_cast<__m256i const*>(b.cell[2 * k + 1]))),
                                       32);
            x[k] = _mm512_castsi512_pd(_mm512_or_si512(lo, hi));
        }
        // cell pairs (int k, int k + 1) as qwords: mid (stack copy for the per-lane old compare)
        alignas(64) std::uint64_t midQ[qC][blockSize];
        __mmask8 post = 0;
        for (std::size_t k = 0; k < qC; ++k)
        {
            __m512i const lo = _mm512_cvtepu32_epi64(
                _mm256_loadu_si256(reinterpret_cast<__m256i const*>(b.midCell[2 * k])));
            __m512i hi = _mm512_setzero_si512();
            if (2 * k + 1 < dim)
                hi = _mm512_slli_epi64(_mm512_cvtepu32_epi64(_mm256_loadu_si256(
                                           reinterpret_cast<__m256i const*>(b.midCell[2 * k + 1]))),
                                       32);
            __m512i const m = _mm512_or_si512(lo, hi);
            _mm512_store_si512(midQ[k], m);
            post |= _mm512_cmpneq_epi64_mask(m, _mm512_castpd_si512(x[k]));
        }
        unsigned pre = 0;
        for (std::size_t l = 0; l < n; ++l)
        {
            auto const* c = reinterpret_cast<std::uint32_t const*>(&particles[first + l].iCell[0]);
            bool diff     = false;
            for (std::size_t k = 0; k < qC; ++k)
            {
                std::uint64_t old = c[2 * k];
                if (2 * k + 1 < dim)
                    old |= std::uint64_t{c[2 * k + 1]} << 32;
                diff |= old != midQ[k][l];
            }
            pre |= unsigned{diff} << l;
        }
        unsigned const live = (1u << n) - 1u;
        for (unsigned m = (pre | post) & live; m; m &= m - 1)
        {
            auto const l = static_cast<std::size_t>(__builtin_ctz(m));
            std::array<int, dim> midCell, newCell;
            for (std::size_t iDim = 0; iDim < dim; ++iDim)
            {
                midCell[iDim] = b.midCell[iDim][l];
                newCell[iDim] = b.cell[iDim][l];
            }
            if (pre >> l & 1u)
                preCrossings_.push_back({first + l, particles[first + l].iCell, midCell});
            if (post >> l & 1u)
                postCrossings_.push_back({first + l, midCell, newCell});
        }
        for (std::size_t iDim = 0; iDim < dim; ++iDim)
            x[qC + iDim] = _mm512_loadu_pd(b.delta[iDim]);
        for (std::size_t c = 0; c < 3; ++c)
            x[qC + dim + c] = _mm512_loadu_pd(b.v[c]);
        for (std::size_t q = nq; q < blockSize; ++q)
            x[q] = _mm512_setzero_pd();
        __m512d t[blockSize], u[blockSize];
        for (std::size_t p = 0; p < blockSize / 2; ++p)
        {
            t[2 * p]     = _mm512_unpacklo_pd(x[2 * p], x[2 * p + 1]);
            t[2 * p + 1] = _mm512_unpackhi_pd(x[2 * p], x[2 * p + 1]);
        }
        for (std::size_t h = 0; h < 2; ++h)
        {
            u[4 * h + 0] = _mm512_shuffle_f64x2(t[4 * h], t[4 * h + 2], 0x88);
            u[4 * h + 1] = _mm512_shuffle_f64x2(t[4 * h], t[4 * h + 2], 0xDD);
            u[4 * h + 2] = _mm512_shuffle_f64x2(t[4 * h + 1], t[4 * h + 3], 0x88);
            u[4 * h + 3] = _mm512_shuffle_f64x2(t[4 * h + 1], t[4 * h + 3], 0xDD);
        }
        auto row = [&](std::size_t const l) __attribute__((always_inline)) {
            constexpr int src[8] = {0, 2, 1, 3, 0, 2, 1, 3};
            return l < 4 ? _mm512_shuffle_f64x2(u[src[l]], u[4 + src[l]], 0x88)
                         : _mm512_shuffle_f64x2(u[src[l]], u[4 + src[l]], 0xDD);
        };
        for (std::size_t l = 0; l < blockSize; ++l)
        {
            if (l >= n)
                break;
            double* const p = reinterpret_cast<double*>(&particles[first + l]) + 1;
            __m512d const y = row(l);
            if constexpr (nq == 8)
                _mm512_storeu_pd(p, y);
            else
            {
                _mm256_storeu_pd(p, _mm512_castpd512_pd256(y));
                _mm_storeu_pd(p + 4, _mm512_extractf64x2_pd(y, 2));
            }
        }
    }

#endif // PHARE_HYB_XSIMD

    /** the particles [begin, end) pushed by blocks of blockSize into block, after(idx, n, ok)
     * called for each block in ascending idx (ok: pushBlock_'s return).
     * By rows on a sorted ParticleArray (CSR valid, whole array) with the cell kernels on, with
     * runs = true (cell deposit of the predictor), or on a bag layout (holes between the rows):
     * blocks per row, the last block of a row partial, then blocks over the unsorted appends.
     * Row blocks take pushBlockCell_ with the cell kernels, else pushBlock_.
     */
    template<typename Interpolator_t, typename Backend, typename Particles, typename After>
    [[gnu::always_inline]] void pushBlocks_(Particles const& particles, std::size_t const begin,
                                            std::size_t const end, BlockFields const& f,
                                            double const coef1, PushBlock& block, After&& after,
                                            bool const runs = false)
    {
        std::size_t idx = begin;
        if constexpr (requires { particles.bagged(); })
        {
#if PHARE_HYB_XSIMD
            bool const cellGather = cellKernels_;
#else
            bool const cellGather = false;
#endif
            assert(!particles.bagged()
                   or (begin == 0 and end == particles.size() and particles.csrValid()));
            if ((cellGather or runs or particles.bagged()) and begin == 0
                and end == particles.size() and particles.csrValid())
            {
#if PHARE_HYB_XSIMD
                CellCache<nbrPointsSupport(Interpolator_t::interp_order)> cache;
#endif
                auto const& csr       = particles.csr();
                int const* const ends = particles.rowEnds();
                auto const& box       = particles.box();
                int const K           = particles.nbrCells();
                std::array<int, dim> cell; // AMR cell of row k
                for (std::size_t i = 0; i < dim; ++i)
                    cell[i] = box.lower[i];
                for (int k = 0; k < K; ++k)
                {
                    auto const e = static_cast<std::size_t>(ends[k]);
#if PHARE_HYB_XSIMD
                    cache.loaded = false;
#endif
                    for (auto i = static_cast<std::size_t>(csr[k]); i < e; i += blockSize)
                    {
                        auto const n = std::min(blockSize, e - i);
#if PHARE_HYB_XSIMD
                        if (cellGather)
                        {
                            after(i, n,
                                  pushBlockCell_<Interpolator_t>(particles, i, n, cell, cache, f,
                                                                 coef1, block));
                            continue;
                        }
#endif
                        after(i, n,
                              pushBlock_<Interpolator_t, Backend>(particles, i, n, f, coef1,
                                                                  block));
                    }
                    for (int d = static_cast<int>(dim) - 1; d >= 0; --d)
                    {
                        if (++cell[d] <= box.upper[d])
                            break;
                        cell[d] = box.lower[d];
                    }
                }
                idx = static_cast<std::size_t>(csr[K]);
            }
        }
        for (; idx < end; idx += blockSize)
        {
            auto const n = std::min(blockSize, end - idx);
            after(idx, n, pushBlock_<Interpolator_t, Backend>(particles, idx, n, f, coef1, block));
        }
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

    bool cellKernels_ = false; // see cellKernels()
};

} // namespace PHARE::core


#endif
