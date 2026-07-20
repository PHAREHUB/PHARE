#ifndef PHARE_CORE_NUMERICS_EULER_HPP
#define PHARE_CORE_NUMERICS_EULER_HPP


#include "core/data/vecfield/vecfield_component.hpp"
#include "core/numerics/finite_volume_euler/finite_volume_euler_per_field.hpp"

namespace PHARE::core
{
template<typename GridLayout>
class FiniteVolumeEuler
{
    constexpr static auto dimension = GridLayout::dimension;


public:
    FiniteVolumeEuler(GridLayout const& layout)
        : layout_{layout}

    {
    }

    template<typename State, typename Fluxes, typename Pred = UpdateAllCells>
    void operator()(State const& state, State& statenew, Fluxes const& fluxes, double const dt,
                    Pred const& shouldUpdate = {}) const
    {
        auto const fve = FiniteVolumeEulerPerField{layout_, dt};

        auto& rhoVxnew = statenew.rhoV(Component::X);
        auto& rhoVynew = statenew.rhoV(Component::Y);
        auto& rhoVznew = statenew.rhoV(Component::Z);

        auto const& rhoVx = state.rhoV(Component::X);
        auto const& rhoVy = state.rhoV(Component::Y);
        auto const& rhoVz = state.rhoV(Component::Z);

        auto const& rhoVx_fx = fluxes.rhoV_fx(Component::X);
        auto const& rhoVy_fx = fluxes.rhoV_fx(Component::Y);
        auto const& rhoVz_fx = fluxes.rhoV_fx(Component::Z);

        if constexpr (dimension == 1)
        {
            fve(state.rho, statenew.rho, shouldUpdate, fluxes.rho_fx);
            fve(rhoVx, rhoVxnew, shouldUpdate, rhoVx_fx);
            fve(rhoVy, rhoVynew, shouldUpdate, rhoVy_fx);
            fve(rhoVz, rhoVznew, shouldUpdate, rhoVz_fx);
            fve(state.Etot, statenew.Etot, shouldUpdate, fluxes.Etot_fx);
        }

        if constexpr (dimension >= 2)
        {
            auto const& rhoVx_fy = fluxes.rhoV_fy(Component::X);
            auto const& rhoVy_fy = fluxes.rhoV_fy(Component::Y);
            auto const& rhoVz_fy = fluxes.rhoV_fy(Component::Z);

            if constexpr (dimension == 2)
            {
                fve(state.rho, statenew.rho, shouldUpdate, fluxes.rho_fx, fluxes.rho_fy);
                fve(rhoVx, rhoVxnew, shouldUpdate, rhoVx_fx, rhoVx_fy);
                fve(rhoVy, rhoVynew, shouldUpdate, rhoVy_fx, rhoVy_fy);
                fve(rhoVz, rhoVznew, shouldUpdate, rhoVz_fx, rhoVz_fy);
                fve(state.Etot, statenew.Etot, shouldUpdate, fluxes.Etot_fx, fluxes.Etot_fy);
            }
            if constexpr (dimension == 3)
            {
                auto const& rhoVx_fz = fluxes.rhoV_fz(Component::X);
                auto const& rhoVy_fz = fluxes.rhoV_fz(Component::Y);
                auto const& rhoVz_fz = fluxes.rhoV_fz(Component::Z);

                fve(state.rho, statenew.rho, shouldUpdate, fluxes.rho_fx, fluxes.rho_fy, fluxes.rho_fz);
                fve(rhoVx, rhoVxnew, shouldUpdate, rhoVx_fx, rhoVx_fy, rhoVx_fz);
                fve(rhoVy, rhoVynew, shouldUpdate, rhoVy_fx, rhoVy_fy, rhoVy_fz);
                fve(rhoVz, rhoVznew, shouldUpdate, rhoVz_fx, rhoVz_fy, rhoVz_fz);
                fve(state.Etot, statenew.Etot, shouldUpdate, fluxes.Etot_fx, fluxes.Etot_fy, fluxes.Etot_fz);
            }
        }
    }

private:
    GridLayout layout_;
};


} // namespace PHARE::core

#endif
