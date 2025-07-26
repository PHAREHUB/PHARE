#ifndef PHARE_CORE_NUMERICS_EULER_HPP
#define PHARE_CORE_NUMERICS_EULER_HPP

#include "initializer/data_provider.hpp"
#include "core/data/grid/gridlayout_utils.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/numerics/finite_volume_euler/finite_volume_euler_per_field.hpp"

namespace PHARE::core
{
template<typename GridLayout>
class FiniteVolumeEuler : public LayoutHolder<GridLayout>
{
    constexpr static auto dimension = GridLayout::dimension;
    using LayoutHolder<GridLayout>::layout_;

public:
    template<typename State, typename Fluxes>
    void operator()(State const& state, State& statenew, Fluxes const& fluxes,
                    double const dt) const
    {
        auto const fve = FiniteVolumeEulerPerField_ref{*layout_, dt};

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
            fve(state.rho, statenew.rho, fluxes.rho_fx);
            fve(rhoVx, rhoVxnew, rhoVx_fx);
            fve(rhoVy, rhoVynew, rhoVy_fx);
            fve(rhoVz, rhoVznew, rhoVz_fx);
            fve(state.Etot, statenew.Etot, fluxes.Etot_fx);
        }

        if constexpr (dimension >= 2)
        {
            auto const& rhoVx_fy = fluxes.rhoV_fy(Component::X);
            auto const& rhoVy_fy = fluxes.rhoV_fy(Component::Y);
            auto const& rhoVz_fy = fluxes.rhoV_fy(Component::Z);

            if constexpr (dimension == 2)
            {
                fve(state.rho, statenew.rho, fluxes.rho_fx, fluxes.rho_fy);
                fve(rhoVx, rhoVxnew, rhoVx_fx, rhoVx_fy);
                fve(rhoVy, rhoVynew, rhoVy_fx, rhoVy_fy);
                fve(rhoVz, rhoVznew, rhoVz_fx, rhoVz_fy);
                fve(state.Etot, statenew.Etot, fluxes.Etot_fx, fluxes.Etot_fy);
            }
            if constexpr (dimension == 3)
            {
                auto const& rhoVx_fz = fluxes.rhoV_fz(Component::X);
                auto const& rhoVy_fz = fluxes.rhoV_fz(Component::Y);
                auto const& rhoVz_fz = fluxes.rhoV_fz(Component::Z);

                fve(state.rho, statenew.rho, fluxes.rho_fx, fluxes.rho_fy, fluxes.rho_fz);
                fve(rhoVx, rhoVxnew, rhoVx_fx, rhoVx_fy, rhoVx_fz);
                fve(rhoVy, rhoVynew, rhoVy_fx, rhoVy_fy, rhoVy_fz);
                fve(rhoVz, rhoVznew, rhoVz_fx, rhoVz_fy, rhoVz_fz);
                fve(state.Etot, statenew.Etot, fluxes.Etot_fx, fluxes.Etot_fy, fluxes.Etot_fz);
            }
        }
    }
};

} // namespace PHARE::core

#endif
