#ifndef PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATING_HPP
#define PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATING_HPP

#include "core/data/field/field_tiles.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include "interpolator.hpp"

namespace PHARE::core
{

template<std::size_t dim, std::size_t interpOrder, bool atomic_ops = false,
         typename Interpolator_t = Interpolator<dim, interpOrder, atomic_ops>>
class Interpolating
{
public:
    template<typename Particles>
    inline void operator()(Particles const& particles, auto& rhoP, auto& rhoC, auto& flux,
                           auto const& layout, double coef = 1.)
    {
        particleToMesh(particles, layout, rhoP, rhoC, flux, coef);
    }

    template<typename Particles>
    void particleToMesh(Particles const& particles, auto const& layout, auto& rhoP, auto& rhoC,
                        auto& flux, double coef = 1.)
        requires(Particles::layout_mode == LayoutMode::AoSPCTS)
    {
        for (std::size_t tidx = 0; tidx < particles().size(); ++tidx)
        {
            auto [rhop, rhoc, F] = tiles_at(tidx, rhoP, rhoC, flux);
            auto const& rho_lay  = tile_layout(rhoP, tidx);
            auto& pctile         = particles()[tidx];
            auto& cps            = pctile();
            for (auto const& bix : cps.local_box())
                for (auto const& p : cps(bix))
                    interp_.particleToMesh(p, rhop, rhoc, F, rho_lay, coef);
        }
    }

    template<typename Particles>
    void particleToMesh(Particles const& particles, auto const& layout, auto& rhoP, auto& rhoC,
                        auto& flux, double coef = 1.)
        requires(Particles::layout_mode == LayoutMode::AoSMapped)
    {
        interp_(particles, rhoP, rhoC, flux, layout, coef);
    }

    Interpolator_t interp_;
};

template<std::size_t dim, std::size_t interpOrder, bool atomic_ops = false>
struct MomentumTensorInterpolating
{
    using Interpolator_t = MomentumTensorInterpolator<dim, interpOrder, atomic_ops>;

    Interpolator_t interp_;

public:
    template<typename Particles_t>
    inline void operator()(Particles_t& particles, auto& momentumTensor, auto const& layout,
                           double mass = 1.)
        requires(Particles_t::layout_mode == LayoutMode::AoSMapped)
    {
        interp_(particles, momentumTensor, layout, mass);
    }

    template<typename Particles_t>
    inline void operator()(Particles_t& particles, auto& momentumTensor, auto const& layout,
                           double mass = 1.)
        requires(Particles_t::layout_mode == LayoutMode::AoSPCTS)
    {
        for (std::size_t tidx = 0; tidx < particles().size(); ++tidx)
        {
            auto mt            = tile_at(momentumTensor, tidx);
            auto const& mt_lay = tile_layout(momentumTensor, tidx);

            auto& pctile = particles()[tidx];
            auto& cps    = pctile();
            for (auto const& bix : cps.local_box())
                interp_(cps(bix), mt, mt_lay, mass);
        }
    }
};

} // namespace PHARE::core

#endif /*PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATING_HPP*/
