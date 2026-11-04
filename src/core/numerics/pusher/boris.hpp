#ifndef PHARE_CORE_PUSHER_BORIS_HPP
#define PHARE_CORE_PUSHER_BORIS_HPP


#include "core/errors.hpp"
#include "core/logger.hpp"
#include "core/utilities/box/box.hpp"


#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <utility>
#include <algorithm>
#include <exception>
#include <functional>


namespace PHARE::core::boris
{


struct MoveTwoCellException : std::exception
{
    MoveTwoCellException(double const d, double const v)
        : delta{d}
        , vel{v}
    {
    }

    double delta, vel;
};


// the particle position after half a time step, as {delta, iCell}, the particle is unchanged
template<typename Particle_t, typename Float, std::size_t dim>
auto static next_position(Particle_t const& p, std::array<Float, dim> const& halfDtOverDl)
{
    std::pair<decltype(p.delta), std::array<int, dim>> next;
    auto& [newDelta, newCell] = next;
    for (std::size_t iDim = 0; iDim < dim; ++iDim)
    {
        double const delta = p.delta[iDim] + static_cast<double>(halfDtOverDl[iDim] * p.v[iDim]);

        if (std::abs(delta) > 2)
            throw MoveTwoCellException{delta, p.v[iDim]};

        auto const iCell = static_cast<int>(std::floor(delta));

        newDelta[iDim] = delta - iCell;
        newCell[iDim]  = iCell + p.iCell[iDim];
    }
    return next;
}

template<typename Particle_t, typename Float, std::size_t dim>
auto static advance(Particle_t& p, std::array<Float, dim> const& halfDtOverDl)
{
    auto const [delta, cell] = next_position(p, halfDtOverDl);
    p.delta                  = delta;
    return cell;
}



/** Accelerate the particles in rangeIn and put the new velocity in rangeOut
 */
template<typename Particle, typename ParticleEB, typename Float>
void accelerate(Particle& p, ParticleEB const& eb, Float const coef1)
{
    static constexpr Float one = 1;
    static constexpr Float two = 2;

    auto& [pE, pB]        = eb;
    auto& [pEx, pEy, pEz] = pE;
    auto& [pBx, pBy, pBz] = pB;

    // We now apply the 3 steps of the BORIS PUSHER
    // 1st half push of the electric field
    Float const velx1 = p.v[0] + coef1 * pEx;
    Float const vely1 = p.v[1] + coef1 * pEy;
    Float const velz1 = p.v[2] + coef1 * pEz;

    // preparing variables for magnetic rotation
    Float const rx = coef1 * pBx;
    Float const ry = coef1 * pBy;
    Float const rz = coef1 * pBz;

    Float const rx2  = rx * rx;
    Float const ry2  = ry * ry;
    Float const rz2  = rz * rz;
    Float const rxry = rx * ry;
    Float const rxrz = rx * rz;
    Float const ryrz = ry * rz;

    Float const invDet = one / (one + rx2 + ry2 + rz2);

    // preparing rotation matrix due to the magnetic field
    // m = invDet*(I + r*r - r x I) - I where x denotes the cross product
    Float const mxx = one + rx2 - ry2 - rz2;
    Float const mxy = two * (rxry + rz);
    Float const mxz = two * (rxrz - ry);

    Float const myx = two * (rxry - rz);
    Float const myy = one + ry2 - rx2 - rz2;
    Float const myz = two * (ryrz + rx);

    Float const mzx = two * (rxrz + ry);
    Float const mzy = two * (ryrz - rx);
    Float const mzz = one + rz2 - rx2 - ry2;

    // magnetic rotation
    Float const velx2 = (mxx * velx1 + mxy * vely1 + mxz * velz1) * invDet;
    Float const vely2 = (myx * velx1 + myy * vely1 + myz * velz1) * invDet;
    Float const velz2 = (mzx * velx1 + mzy * vely1 + mzz * velz1) * invDet;

    // 2nd half push of the electric field / Update particle velocity
    p.v[0] = velx2 + coef1 * pEx;
    p.v[1] = vely2 + coef1 * pEy;
    p.v[2] = velz2 + coef1 * pEz;
}


} // namespace PHARE::core::boris




#endif
