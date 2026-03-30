#ifndef PHARE_CORE_PUSHER_BORIS_HPP
#define PHARE_CORE_PUSHER_BORIS_HPP


#include "core/errors.hpp"
#include "core/logger.hpp"
#include "core/utilities/box/box.hpp"


#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
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


template<typename Particle_t, typename Float, std::size_t dim>
auto static advance(Particle_t& p, std::array<Float, dim> const& halfDtOverDl)
{
    std::array<int, dim> newCell;
    for (std::size_t iDim = 0; iDim < dim; ++iDim)
    {
        double const delta = p.delta[iDim] + static_cast<double>(halfDtOverDl[iDim] * p.v[iDim]);

        if (std::abs(delta) > 2)
            throw MoveTwoCellException{delta, p.v[iDim]};

        auto const iCell = static_cast<int>(std::floor(delta));

        p.delta[iDim] = delta - iCell;
        newCell[iDim] = iCell + p.iCell[iDim];
    }
    return newCell;
}



/** Accelerate the particles in rangeIn and put the new velocity in rangeOut
 */
template<typename Particle, typename ParticleEB, typename Float>
void accelerate(Particle& p, ParticleEB const& eb, Float const& dto2m)
{
    static constexpr Float one = 1;
    static constexpr Float two = 2;

    auto& [pE, pB]        = eb;
    auto& [pEx, pEy, pEz] = pE;
    auto& [pBx, pBy, pBz] = pB;

    Float const coef1 = p.charge * dto2m;

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


namespace PHARE::core
{


template<std::size_t dim, typename ParticleArray_t, typename Electromag, typename Interpolator,
         typename BoundaryCondition, typename GridLayout>
class BorisPusher
{
    using Particle_t    = ParticleArray_t::value_type;
    using OnChangeIcell = std::function<bool(ParticleArray_t&, Particle_t&,
                                             std::array<int, dim> const&, std::size_t const)>;


public:
    void setMeshAndTimeStep(std::array<double, dim> ms, double const ts)
    {
        std::transform(std::begin(ms), std::end(ms), std::begin(halfDtOverDl_),
                       [ts](double& x) { return 0.5 * ts / x; });
        dt_ = ts;
    }


    template<bool copy_particle = false, typename Population, typename Boxing_t>
    void move(Population& pop, Electromag const& em, Interpolator& interpolator,
              Boxing_t const& boxing)
    {
        move_domain<copy_particle>(pop, em, interpolator, boxing);
        move_level_ghost<copy_particle>(pop, em, interpolator, boxing);
    }

private:
    void move_particle(auto& particle, auto const& em, auto& interpolator, auto const& layout,
                       auto const& dto2m)
    {
        try
        {
            particle.iCell = boris::advance(particle, halfDtOverDl_);
        }
        catch (boris::MoveTwoCellException const& e)
        {
            std::stringstream ss;
            ss << "PrePush Particle moved 2 cells with delta/vel: ";
            ss << e.delta << "/" << e.vel << std::endl;
            throw DictionaryException{ss.str()};
        }

        auto const& local_em = interpolator(particle, em, layout);
        try
        {
            boris::accelerate(particle, local_em, dto2m);
            particle.iCell = boris::advance(particle, halfDtOverDl_);
        }
        catch (boris::MoveTwoCellException const& e)
        {
            std::stringstream ss;
            ss << "PostPush Particle moved 2 cells with delta/vel: ";
            ss << e.delta << "/" << e.vel << std::endl;
            DictionaryException ex{ss.str()};

            auto const& [E, B] = local_em;
            for (std::uint16_t i = 0; i < 3; ++i)
                ex("E_" + std::to_string(i), std::to_string(E[i]));
            for (std::uint16_t i = 0; i < 3; ++i)
                ex("B_" + std::to_string(i), std::to_string(B[i]));
            ex("level", std::to_string(layout.levelNumber()));
            throw ex;
        }
    };

    template<bool copy_particle, typename Population, typename Boxing_t>
    void move_level_ghost(Population& pop, Electromag const& em, Interpolator& interpolator,
                          Boxing_t const& boxing)
    {
        PHARE_LOG_SCOPE(3, "Boris::move_level_ghost");

        using ParticleLoop_t = std::conditional_t<copy_particle, Particle_t, Particle_t&>;

        auto const& layout = boxing.layout;
        double const dto2m = 0.5 * dt_ / pop.mass();

        auto& domain      = pop.domainParticles();
        auto& level_ghost = pop.levelGhostParticles();

        for (std::size_t i = 0; i < level_ghost.size(); ++i) // size might change on iteration!
        {
            ParticleLoop_t particle = level_ghost[i];
            auto const oldCell      = particle.iCell;

            move_particle(particle, em, interpolator, layout, dto2m);

            if (oldCell == particle.iCell)
                continue;

            bool const isInDomainBox      = boxing.isInDomainBox(particle);
            bool const should_interpolate = isInDomainBox;

            if (should_interpolate)
                interpolator.particleToMesh( //
                    particle, pop.particleDensity(), pop.chargeDensity(), pop.flux(), layout);

            if constexpr (!copy_particle)
            {
                if (isInDomainBox)
                {
                    domain.push_back(particle);
                    level_ghost.swap_last_reduce_by_one(oldCell, i);
                    --i; // redo current index as last is now i
                    continue;
                }

                bool const isInNonLevelGhostBox
                    = isInDomainBox || boxing.isInNonLevelGhostBox(particle.iCell);
                bool const isInLevelGhostBox = !isInNonLevelGhostBox;

                if (isInLevelGhostBox)
                    level_ghost.change_icell(particle, oldCell, i);
                else
                {
                    level_ghost.swap_last_reduce_by_one(oldCell, i);
                    --i; // redo current index as last is now i
                }
            }
        }
    }

    template<bool copy_particle, typename Population, typename Boxing_t>
    void move_domain(Population& pop, Electromag const& em, Interpolator& interpolator,
                     Boxing_t const& boxing)
    {
        PHARE_LOG_SCOPE(3, "Boris::move_domain");

        using ParticleLoop_t = std::conditional_t<copy_particle, Particle_t, Particle_t&>;


        auto const& layout = boxing.layout;
        double const dto2m = 0.5 * dt_ / pop.mass();

        auto& domain      = pop.domainParticles();
        auto& patch_ghost = pop.patchGhostParticles();

        for (std::size_t i = 0; i < domain.size(); ++i) // size might change on iteration!
        {
            ParticleLoop_t particle = domain[i];
            auto const oldCell      = particle.iCell;

            move_particle(particle, em, interpolator, layout, dto2m);

            if (oldCell == particle.iCell)
            {
                interpolator.particleToMesh( //
                    particle, pop.particleDensity(), pop.chargeDensity(), pop.flux(), layout);

                continue;
            }

            bool const isInDomainBox = boxing.isInDomainBox(particle);
            bool const isInNonLevelGhostBox
                = isInDomainBox || boxing.isInNonLevelGhostBox(particle.iCell);
            bool const should_interpolate = isInNonLevelGhostBox;

            if (!should_interpolate)
            {
                if constexpr (!copy_particle)
                {
                    domain.swap_last_reduce_by_one(oldCell, i);
                    --i; // redo current index as last is now i
                }
                continue;
            }

            interpolator.particleToMesh( //
                particle, pop.particleDensity(), pop.chargeDensity(), pop.flux(), layout);

            if constexpr (!copy_particle)
                domain.change_icell(particle, oldCell, i);
        }

        // move to patch_ghost
        if constexpr (!copy_particle)
        {
            auto range = makeIndexRange(domain);
            range      = domain.partition(
                range, [&](auto const& cell) { return core::isIn(cell, boxing.domainBox); });

            auto const not_in_domain  = makeRange(domain, range.iend(), domain.size());
            auto const is_patch_ghost = domain.partition(
                not_in_domain, [&](auto const& cell) { return boxing.isInNonLevelGhostBox(cell); });

            patch_ghost.reserve(patch_ghost.size() + is_patch_ghost.size());
            std::copy(is_patch_ghost.begin(), is_patch_ghost.end(),
                      std::back_inserter(patch_ghost));
            domain.erase(not_in_domain);
        }

        PHARE_DEBUG_DO({
            for (auto const& p : domain)
                if (!isIn(p, boxing.domainBox))
                    throw std::runtime_error("invalid domain");
            for (auto const& p : patch_ghost)
                if (isIn(p, boxing.domainBox))
                    throw std::runtime_error("invalid patch ghost");
        })
    }

    std::array<double, dim> halfDtOverDl_;
    double dt_;
};



} // namespace PHARE::core


#endif
