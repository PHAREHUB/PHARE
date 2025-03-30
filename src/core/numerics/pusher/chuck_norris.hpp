#ifndef PHARE_CORE_PUSHER_CHUCK_NORRIS_HPP
#define PHARE_CORE_PUSHER_CHUCK_NORRIS_HPP


#include "core/errors.hpp"
#include "core/logger.hpp"

#include "core/numerics/pusher/fusher.hpp"
#include "core/data/particles/particle.hpp"


#include <array>
#include <cmath>
#include <core/utilities/types.hpp>
#include <cstddef>
#include <algorithm>
#include <iterator>
#include <stdexcept>

namespace PHARE::core::other
{


template<std::size_t dim, typename ParticleArray_t, typename Electromag, typename Interpolator,
         typename BoundaryCondition, typename GridLayout>
class BorisPusher
    : public Pusher<dim, ParticleArray_t, Electromag, Interpolator, BoundaryCondition, GridLayout>
{
public:
    using Super
        = Pusher<dim, ParticleArray_t, Electromag, Interpolator, BoundaryCondition, GridLayout>;

private:
    using OnChangeIcell = typename Super::OnChangeIcell;

public:
    void move(ParticleArray_t& particles, Electromag const& emFields, double mass,
              Interpolator& interpolator, GridLayout const& layout, OnChangeIcell&& OnChangeIcell,
              bool const ghost) override
    {
        PHARE_LOG_SCOPE(3, "Boris::move_no_bc");
        auto constexpr partGhostWidth = GridLayout::nbrParticleGhosts();

        [[maybe_unused]] auto const ghostBox = grow(layout.AMRBox(), partGhostWidth);
        double const dto2m                   = 0.5 * dt_ / mass;

        auto const accer = [&](auto& particle) {
            PHARE_LOG_LINE_SS(Point{particle.iCell});
            accelerate_(particle, interpolator(particle, emFields, layout), dto2m);
            particle.iCell = advancePosition_(particle);
        };

        auto const mover = [&]<bool GHOST>() {
            for (std::size_t i = 0; i < particles.size(); ++i) // size might change on iteration!
            {
                auto& particle = particles[i];
                PHARE_LOG_LINE_SS(Point{particle.iCell});
                auto const oldCell = particle.iCell;
                particle.iCell     = advancePosition_(particle);

                if constexpr (GHOST)
                {
                    if (isIn(particle, ghostBox))
                        accer(particle);
                }
                else
                    accer(particle);

                if (oldCell != particle.iCell)
                    if (OnChangeIcell(particles, particle, oldCell, i))
                        --i; // repeat current index because of swap last to idx
            }
        };

        if (ghost)
            mover.template operator()<1>();
        else
            mover.template operator()<0>();
    }



    /** see Pusher::move() documentation*/
    void setMeshAndTimeStep(std::array<double, dim> ms, double const ts) override
    {
        std::transform(std::begin(ms), std::end(ms), std::begin(halfDtOverDl_),
                       [ts](double& x) { return 0.5 * ts / x; });
        dt_ = ts;
    }



private:
    /** move the particle partIn of half a time step and store it in partOut
     */
    template<typename Particle>
    auto advancePosition_(Particle& particle)
    {
        std::array<int, dim> newCell;
        for (std::size_t iDim = 0; iDim < dim; ++iDim)
        {
            double const delta = particle.delta[iDim]
                                 + static_cast<double>(halfDtOverDl_[iDim] * particle.v[iDim]);

            double const iCell = std::floor(delta);
            if (std::abs(delta) > 2)
            {
                PHARE_LOG_ERROR("Error, particle moves more than 1 cell, delta >2");
            }
            particle.delta[iDim] = delta - iCell;
            newCell[iDim]        = static_cast<int>(iCell + particle.iCell[iDim]);
        }
        return newCell;
    }




    /** Accelerate the particles in rangeIn and put the new velocity in rangeOut
     */
    template<typename Particle_t, typename ParticleEB>
    void accelerate_(Particle_t& part, ParticleEB const& particleEB, double const& dto2m)
    {
        auto& [pE, pB]        = particleEB;
        auto& [pEx, pEy, pEz] = pE;
        auto& [pBx, pBy, pBz] = pB;


        double const coef1 = part.charge * dto2m;

        // We now apply the 3 steps of the BORIS PUSHER

        // 1st half push of the electric field
        double const velx1 = part.v[0] + coef1 * pEx;
        double const vely1 = part.v[1] + coef1 * pEy;
        double const velz1 = part.v[2] + coef1 * pEz;


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


        // 2nd half push of the electric field / Update particle velocity
        part.v[0] = velx2 + coef1 * pEx;
        part.v[1] = vely2 + coef1 * pEy;
        part.v[2] = velz2 + coef1 * pEz;
    }




    std::array<double, dim> halfDtOverDl_;
    double dt_;
};

} // namespace PHARE::core::other


#endif
