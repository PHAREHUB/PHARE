#ifndef PHARE_CORE_PUSHER_BORIS_HPP
#define PHARE_CORE_PUSHER_BORIS_HPP


#include "core/errors.hpp"
#include "core/logger.hpp"
#include "core/numerics/pusher/pusher.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/models/quantities/hybrid_quantities.hpp"


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
        for (auto idx = range.ibegin(); idx < end; idx += blockSize)
        {
            auto const n = std::min(blockSize, end - idx);
            if (!pushBlock_<Interpolator>(particles, idx, n, fields, coef1, block)) [[unlikely]]
            {
                // a particle of the block moves more than one cell: the scalar path throws
                // the same exception as move()
                for (auto i = idx; i < idx + n; ++i)
                    moveOne(i);
                continue;
            }
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


    /** push() on a copy of every particle of particles, in order, calling onPushed(pushed)
     * after each. Particles are pushed by blocks of blockSize with the same per-particle
     * arithmetic as push(); a block where a particle moves more than one cell is redone with
     * push(), which throws.
     */
    template<typename Particles, typename OnPushed>
    [[gnu::noinline]] void pushEach(Particles const& particles, Electromag const& emFields,
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
        for (std::size_t idx = 0; idx < end; idx += blockSize)
        {
            auto const n = std::min(blockSize, end - idx);
            if (!pushBlock_<Interpolator>(particles, idx, n, fields, coef1, block)) [[unlikely]]
            {
                for (auto i = idx; i < idx + n; ++i)
                    pushOne(i);
                continue;
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
        }
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

    static void accelerate_(double& vx, double& vy, double& vz, double const pEx,
                            double const pEy, double const pEz, double const pBx,
                            double const pBy, double const pBz, double const coef1)
    {
        // We now apply the 3 steps of the BORIS PUSHER

        // 1st half push of the electric field
        double velx1 = vx + coef1 * pEx;
        double vely1 = vy + coef1 * pEy;
        double velz1 = vz + coef1 * pEz;


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
    static bool advanceBlock_(std::array<double, dim> const& halfDtOverDl,
                              double (&delta)[dim][blockSize], int (&cell)[dim][blockSize],
                              double const (&v)[3][blockSize])
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
    static void gatherBlock_(BlockField const& f, Starts const& starts, Weights const& weights,
                             double (&out)[blockSize])
    {
        double const* __restrict__ data = f.data;
        if constexpr (dim == 1)
        {
            auto const& [xs]    = starts;
            auto const& [xw]    = weights;
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
    template<typename Interpolator_t, QtyCentering centering, std::size_t nbrPts>
    static void indexAndWeightsBlock_(int const (&localCell)[dim][blockSize],
                                      double const (&delta)[dim][blockSize],
                                      int (&starts)[dim][blockSize],
                                      double (&weights)[dim][nbrPts][blockSize])
    {
        constexpr auto order = Interpolator_t::interp_order;
        for (std::size_t iDim = 0; iDim < dim; ++iDim)
            for (std::size_t l = 0; l < blockSize; ++l)
            {
                double const d = delta[iDim][l];
                int const lc   = localCell[iDim][l];
                int const s    = lc
                              - Interpolator_t::template computeStartLeftShift<QtyCentering,
                                                                               centering>(d);
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
    template<typename Interpolator_t, typename Particles>
    bool pushBlock_(Particles const& particles, std::size_t const first, std::size_t const n,
                    BlockFields const& f, double const coef1, PushBlock& b) const
    {
        constexpr std::size_t nbrPts = nbrPointsSupport(Interpolator_t::interp_order);
        using Scalar                 = HybridQuantity::Scalar;

        auto const halfDtOverDl = halfDtOverDl_;

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

        bool tooFar = advanceBlock_(halfDtOverDl, b.delta, b.cell, b.v);

        int localCell[dim][blockSize];
        for (std::size_t iDim = 0; iDim < dim; ++iDim)
            for (std::size_t l = 0; l < blockSize; ++l)
            {
                b.midCell[iDim][l]   = b.cell[iDim][l];
                localCell[iDim][l] = b.cell[iDim][l] + f.amrToLocal[iDim];
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

        for (std::size_t l = 0; l < blockSize; ++l)
            accelerate_(b.v[0][l], b.v[1][l], b.v[2][l], EB[0][l], EB[1][l], EB[2][l], EB[3][l],
                        EB[4][l], EB[5][l], coef1);

        tooFar |= advanceBlock_(halfDtOverDl, b.delta, b.cell, b.v);
        return !tooFar;
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
